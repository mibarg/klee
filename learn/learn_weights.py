
import subprocess
import argparse
from os import walk
from os.path import dirname, join, abspath
import json
from random import random
import pandas as pd
from datetime import datetime
import builtins
from scipy.optimize import basinhopping, minimize
from stats import main as get_stats
import numpy as np

EPS = 0.0001
CLR = '\033[1;34m'  # Light Blue, see: https://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
NCLR = '\033[0m' # No Color
LOGFILE = None
F_CACHE = {}

features = [
    "instr_pc_getNumOperands", "instr_pc_getOpcode", "instr_pc_use_empty", "instr_prevpc_getNumOperands",
    "instr_prevpc_getOpcode", "instr_prevpc_use_empty", "sat_paper_num_ands", "sat_paper_num_asserted",
    "sat_paper_num_eq", "sat_paper_num_xor", "sat_paper_num_bool", "sat_paper_num_theory", "sat_paper_clause_to_vars",
    "sat_paper_vars_clause", "sat_org_num_sub", "sat_org_num_mul", "sat_org_num_udiv", "sat_org_num_sdiv",
    "sat_org_num_urem", "sat_org_num_or", "sat_org_num_shl", "sat_org_num_lshr", "sat_org_num_ashr",
    "sat_org_num_extract", "sat_org_num_concat", "sat_org_num_zext", "sat_org_num_sext", "sat_org_num_lsb",
    "sat_org_num_msb", "sat_org_cls_std", "sat_org_cls_avg", "sat_org_cls_max", "sat_org_num_capitals", 
    "weights_Depth", "weights_InstCount", "weights_CPInstCount", "weights_QueryCost", "weights_CoveringNew",
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
    
    # log to file
    if LOGFILE is not None:
        LOGFILE.write('%s\n' % ' '.join(nargs))
    
    return builtins.print(*nargs, **kwargs)

    
def parse_args(argv=None):
    """
    Define and parse scipt arguments
    :param argv: if provided, these will be used instead of sys.argv
    """

    parser = argparse.ArgumentParser()

    # klee
    parser.add_argument('-kd', "--klee-dir", default='', help="KLEE directory")
    parser.add_argument('-ks', "--klee-searcher", default='smart-weighted', help="KLEE searcher name")
    parser.add_argument('-ki', "--klee-input", default='', help="KLEE input file")
    parser.add_argument('-ka', "--klee-args", default='', help="KLEE args")
    parser.add_argument('-kb', "--klee-args-after", default='', help="KLEE args that need to come after the input file (ki).")

    # goal
    parser.add_argument('-g', "--goal", default='CoveredInstructions', help="KLEE stat or NumBugs to maximize")
    parser.add_argument('-w', "--weights", default='weights.json', help="Weights for KLEE searcher, in JSON format")
    parser.add_argument('-nr', "--num-runs", default=2, type=int, help="Number of KLEE sessions to run")
    parser.add_argument('-ni', "--num-inits", default=100, type=int, help="Number of random initial values x0 to try for each optimization session")
    parser.add_argument('-wi', "--weights-init", default=False, action='store_true', help="If True, initiate weights vector to contain only weight_* features")
    parser.add_argument('-si', "--save-interval", default=1, type=int, help="Save data interval (of KLEE sessions)")
    
    # output
    parser.add_argument('-on', "--output-name", default='oracle_cache.pckl', help="Output filename")
    parser.add_argument('-s', "--silent", default=False, action='store_true', help="No KLEE output logs")
    parser.add_argument('-ln', "--log-name", default='log.txt', help="Log file")
    
    # scipy
    parser.add_argument('-mm', "--minimizer-method", default=None, help="Kwargs for scipy minimizer")
    
    if argv is None:
        return parser.parse_args()
    else:
        return parser.parse_args(argv)


def run_klee_session(args):
    """
    Run a KLEE session
    """
    # import pdb
    # pdb.set_trace()
    cmnd = '%s/klee --search=%s %s %s %s' % (args.klee_dir, args.klee_searcher,   args.klee_args, args.klee_input, args.klee_args_after)

    try:
        if args.silent:
            #print('Running silent klee session with cmnd: %s' % cmnd)
            proc = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            #print('Running noisy klee session with cmnd: %s' % cmnd)
            proc = subprocess.Popen(cmnd, shell=True)
        # FIXME choose
        proc.communicate() # or wait
        return proc.returncode
    except Exception:
        return -1


def read_stat(args):
    """
    Read args.goal from stats file
    Count number of .err files if args.goal == 'NumBugs'
    Read coverage % if args.goal in ['BCov(%)', 'ICov(%)']
    """
    klee_last_path = join(dirname(abspath(args.klee_input)), 'klee-last')
    if args.goal == 'NumBugs':
        # count .err files
        return float(len([file for curdir, subdirs, file in walk(klee_last_path) if '.err' in file]))
    if args.goal in ['BCov(%)', 'ICov(%)']:
        stats = get_stats([klee_last_path,])
        return float(stats[1][stats[0].index(args.goal)])
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

                return float(values[cols.index(args.goal)])

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
		

def set_weights(weights, args):
    """
    Write weights file
    """
    
    weights_path = join(dirname(args.klee_input), args.weights)
    try:
        # save weights
        with open(weights_path, 'w') as f:
            f.write(json.dumps(weights))

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

    # random weights
    new_weights = {feat: random() * 2 - 1 for feat in features}

    # save new weights
    return set_weights(new_weights, args)
		

def func(w, args):
    global F_CACHE

    # check cache
    if tuple(w) in F_CACHE:
        return F_CACHE[tuple(w)]

    # set weights JSON
    weights_dict = {features[i]: w[i] for i in range(len(features))}
    if set_weights(weights_dict, args) != 0: return None

    # run searcher
    if run_klee_session(args) != 0:  return None

    # read cov % from stats
    goal_value = read_stat(args)
    if goal_value == -1:  return None
    
    # update cache
    F_CACHE[tuple(w)] = goal_value
    
    print('f(w) = %.4f' % goal_value)
    return goal_value
    
    
def learn_by_grad(args):
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
    
    
def define_logger(args):
    """
    define logfile, if requested by args
    """
    global LOGFILE
    logfile_path = join(dirname(args.klee_input), args.log_name)
    LOGFILE = open(logfile_path, 'w')
    print('Saving logfile to %s' % abspath(logfile_path))

    
def main():
    global F_CACHE

    args = parse_args()
    define_logger(args)
    
    print('########## learn_weights start ##########')
    print('Preparing to run %d KLEE sessions' % args.num_runs)
    
    # goal function
    f = lambda w: -func(w, args)
    if args.weights_init:
        x0 = [1.0 if 'weights_' in feat else EPS for feat in features]
    else:
        x0 = [random() * 2 - 1 for feat in features]
    callb = lambda x: print("minimum %.4f, found at %s" % (None if x is None else f(x), x))
    
    try:
        if args.minimizer_method == 'random':
            best_w, best_v = None, 0
            for i in range(args.num_runs):
                rnd_w = [random() * 2 - 1 for feat in features]
                cur_v = f(rnd_w)
                if cur_v < best_v:
                    best_v = cur_v
                    best_w = rnd_w
            print('Global minimum found: %.4f' % best_v)
            ret_w = {features[i]: best_w for i in range(len(features))}
            set_weights(ret_w, args)
        else:
        
            best_w, best_res = None, 0
            # 
            for i in range(args.num_inits):
                #ret = basinhopping(f, x0, 
                #    niter=args.num_runs, 
                #    minimizer_kwargs={'method': args.minimizer_method}, 
                #    callback=callb,
                #    interval=1,
                #)
                #'initial_simplex': np.vstack([np.eye(len(features)) - EPS * np.ones((len(features), len(features))), -1 * np.ones((1,len(features)))]),
                
                # Nelder-Mead, Powell, COBYLA print issue
                tol=1e-7
                ret = minimize(f, x0, 
                    method = args.minimizer_method,
                    callback=callb,
                    tol=tol,
                    options={
                        'disp': True, 
                        'maxiter': args.num_runs, 
                        'xtol': tol, 'ftol': tol, 'xatol': tol, 'fatol': tol}
                )
                
                if ret.fun < best_res:
                    best_w = ret.x
                    best_res = ret.fun
                
            # report and save result
            print('Global minimum found: %.4f' % best_res)
            
            ret_w = {features[i]: best_w[i] for i in range(len(features))}
            set_weights(ret_w, args)
        
            # print configuration
            print(ret)
    except KeyboardInterrupt:
        print('Optimization interrupted')
    
    # save F_CACHE
    pd.DataFrame \
        .from_dict(F_CACHE, orient='index') \
        .reset_index().set_index(0) \
        ['index'].apply(pd.Series) \
        .to_pickle(args.output_name)
    print('F_CACHE saved to %s' % abspath(args.output_name))
    
    print('########## learn_weights end ##########')
    if hasattr(LOGFILE, 'close'): LOGFILE.close()


if __name__ == '__main__':
    main()