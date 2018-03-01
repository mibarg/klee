
import subprocess
import argparse
from os import walk
from os.path import dirname, join
import json
from random import random
import pandas as pd
from datetime import datetime
import builtins

EPS = 0.0001
CLR = '\033[1;34m'  # Light Blue, see: https://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
NCLR = '\033[0m' # No Color

features = [
    "instr_pc_getNumOperands", "instr_pc_getOpcode", "instr_pc_use_empty", "instr_prevpc_getNumOperands",
    "instr_prevpc_getOpcode", "instr_prevpc_use_empty", "sat_paper_num_ands", "sat_paper_num_asserted",
    "sat_paper_num_eq", "sat_paper_num_xor", "sat_paper_num_bool", "sat_paper_num_theory", "sat_paper_clause_to_vars",
    "sat_paper_vars_clause", "sat_org_num_sub", "sat_org_num_mul", "sat_org_num_udiv", "sat_org_num_sdiv",
    "sat_org_num_urem", "sat_org_num_or", "sat_org_num_shl", "sat_org_num_lshr", "sat_org_num_ashr",
    "sat_org_num_extract", "sat_org_num_concat", "sat_org_num_zext", "sat_org_num_sext", "sat_org_num_lsb",
    "sat_org_num_msb", "sat_org_cls_std", "sat_org_cls_avg", "sat_org_cls_max", "sat_org_num_capitals", "weights_Depth",
    "weights_InstCount", "weights_CPInstCount", "weights_QueryCost", "weights_CoveringNew",
    "weights_MinDistToUncovered", "weights_forkDisabled"]

    

def print(*args, **kwargs):
    """
    Customize prints
    """
    frst_arg = '%s[%s] LEARN: %s%s' % (
            CLR, 
            datetime.now().strftime("%Y/%m/%d %H:%M:%S"), 
            args[0], NCLR)
    if len(args) > 1:
        nargs = tuple([frst_arg, list(args[1:])])
    else:
        nargs = (frst_arg,)
    return builtins.print(*nargs, **kwargs)

    
def parse_args(argv=None):
    """
    Define and parse scipt arguments
    :param argv: if provided, these will be used instead of sys.argv
    """

    parser = argparse.ArgumentParser()

    # klee
    parser.add_argument('-kd', "--klee-dir", default='', help="KLEE directory")
    parser.add_argument('-ks', "--klee-searcher", default='smart_weighted', help="KLEE searcher name")
    parser.add_argument('-ki', "--klee-input", default='', help="KLEE input file")
    parser.add_argument('-ka', "--klee-args", default='', help="KLEE args")

    # goal
    parser.add_argument('-g', "--goal", choices=('CoveredInstructions', 'NumBugs'), default='CoveredInstructions', help="KLEE stat or NumBugs to maximize")
    parser.add_argument('-w', "--weights", default='weights.json', help="Weights for KLEE searcher, in JSON format")
    parser.add_argument('-nr', "--num-runs", default=2, type=int, help="Number of KLEE sessions to run")
    parser.add_argument('-si', "--save-interval", default=1, type=int, help="Save data interval (of KLEE sessions)")
    
    parser.add_argument('-on', "--output-name", default='weights_to_goal.pckl', help="Output filename")
    parser.add_argument('-s', "--silent", default=False, action='store_true', help="No KLEE output logs")

    # examples
    #parser.add_argument('-l', "--logging-level", choices=('NOTSET', 'INFO', 'DEBUG', 'ERROR'), default='INFO', help="Level of logging messages to show")
    #parser.add_argument('-sp', "--port", default=7070, type=int, help="Socket server port")
    #parser.add_argument('-r', "--restore", default=False, action='store_true', help="Restore network parameters from last save, if they exist")

    if argv is None:
        return parser.parse_args()
    else:
        return parser.parse_args(argv)


def run_klee_session(args):
    """
    Run a KLEE session
    """

    cmnd = '%s/klee --search=%s %s %s' % (args.klee_dir, args.klee_searcher, args.klee_args, args.klee_input)

    try:
        if args.silent:
            print('Running silent klee session with cmnd: %s' % cmnd)
            proc = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE)
        else:
            print('Running noisy klee session with cmnd: %s' % cmnd)
            proc = subprocess.Popen(cmnd, shell=True)
        # FIXME choose
        proc.communicate() # or wait
        return proc.returncode
    except Exception:
        return -1


def read_stat(args):
    """
    Read args.goal from stats file, or count number of .err files if args.goal == 'NumBugs'
    """

    if args.goal == 'NumBugs':
        # count .err files
        return float(len([file for curdir, subdirs, file in walk('klee-last/') if '.err' in file]))
    else:
        # search .stat file
        stats_path = join(dirname(args.klee_input), 'klee-last', 'run.stats')
        try:
            # get relative dir
            with open(stats_path, 'r') as f:
                # parse header
                header = f.readline().strip()
                cols = header.replace("('", "").replace("',)", "").split("','")
                goal_ind = cols.index(args.goal)

                # iterate to last line
                line = ''
                for line in f: pass
                clean_line = line.strip().replace("(", "").replace(")", "")
                values = clean_line.split(",")
				
                covered_inst = float(values[cols.index('CoveredInstructions')])
                uncovered_inst = float(values[cols.index('UncoveredInstructions')])

                return covered_inst / (covered_inst + uncovered_inst + EPS)

        except (ValueError, IndexError, OSError) as e:
            if isinstance(e, (ValueError, IndexError)):
                print('Error while looking for %s in %s: %s' % (args.goal, stats_path, e))
            if isinstance(e, OSError):
                print('Error on read %s: %s' % (stats_path, e))
            return -1


def gradient_step(weights, value):
    return weights


def update_weights(args, goal_value, type='gradient'):
    """
    Update weights file based on goal_value
    type: can be 'gradient' or 'random'
    """
    assert(type in ('gradient', 'random'))

    weights_path = join(dirname(args.klee_input), args.weights)
    try:
        # get current weights
        current_weights = get_weights(args)
        if current_weights is None:
            return -1
        
        if type == 'gradient':
            # gradient step
            new_weights = gradient_step(current_weights, goal_value)

            # save new weights
            with open(weights_path, 'w') as f:
                f.write(json.dumps(new_weights))
        else:
            # type == 'random'
            init_weights(args)

        return 0

    except (ValueError , OSError) as e:
        if isinstance(e, ValueError):
            print('Error while loading/dumping weights to JSON: %s' % e)
        if isinstance(e, OSError):
            print('Error on read/write %s: %s' % (weights_path, e))
        return -1

		
def get_weights(args):
    """
    Read weights file
    """

    weights_path = join(dirname(args.klee_input), args.weights)
    try:
        # get current weights
        current_weights = None
        with open(weights_path, 'r') as f:
            current_weights = json.loads(f.read())

        return current_weights

    except (ValueError , OSError) as e:
        if isinstance(e, ValueError):
            print('Error while loading/dumping weights to JSON: %s' % e)
        if isinstance(e, OSError):
            print('Error on read/write %s: %s' % (weights_path, e))
        return None


def init_weights(args):
    """
    Randomly init a weights file
    """

    weights_path = join(dirname(args.klee_input), args.weights)
    try:
        # random weights
        new_weights = {feat: random() * 2 - 1 for feat in features}

        # save new weights
        with open(weights_path, 'w') as f:
            f.write(json.dumps(new_weights))

        return 0

    except (json.JSONDecodeError, OSError) as e:
        if isinstance(e, json.JSONDecodeError):
            print('Error while loading/dumping weights to JSON: %s' % e)
        if isinstance(e, OSError):
            print('Error on read/write %s: %s' % ( weights_path, e))
        return -1


def main():
    args = parse_args()
    print('########## learn_weights start ##########')
    print('Preparing to run %d KLEE sessions' % args.num_runs)

    failures = dict(klee=0, stats=0, update=0)
    cntr = 0
    data = None

    # random start
    if init_weights(args) != 0:
        print('Error on weight init')
        exit(0)

    # weight iteration
    while True:
        print('KLEE session number %d' % cntr)

        # run searcher
        if run_klee_session(args) != 0:
            # something went wrong with the klee sessions
            failures['klee'] += 1
            continue

        # read cov % from stats
        goal_value = read_stat(args)
        if goal_value == -1:
            # something went wrong with reading run.stats
            failures['stats'] += 1
            continue
			
        # record data
        data_dict = {**get_weights(args), **{'goal': goal_value}}
        if data is None:
            data = pd.DataFrame.from_dict({0: data_dict}, orient='index')
        else:
            data = data.append(data_dict, ignore_index=True)

        # update weights
        if update_weights(args, goal_value, 'random') != 0:
            # something went wrong with the weights update
            failures['update'] += 1
            continue
        
        # save data
        if cntr % args.save_interval == 0:
            data.to_pickle(args.output_name)

        # break condition
        if cntr >= (args.num_runs - 1):
            data.to_pickle(args.output_name)
            break
        else:
            cntr += 1
                
    print('Data collected: (%d, %d)' % (data.shape[0], data.shape[1]))
    print('########## learn_weights end ##########')


if __name__ == '__main__':
    main()

