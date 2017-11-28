//===-- Searcher.cpp ------------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "Searcher.h"

#include "CoreStats.h"
#include "Executor.h"
#include "PTree.h"
#include "StatsTracker.h"

#include "klee/ExecutionState.h"
#include "klee/Statistics.h"
#include "klee/Internal/Module/InstructionInfoTable.h"
#include "klee/Internal/Module/KInstruction.h"
#include "klee/Internal/Module/KModule.h"
#include "klee/Internal/ADT/DiscretePDF.h"
#include "klee/Internal/ADT/RNG.h"
#include "klee/Internal/Support/ModuleUtil.h"
#include "klee/Internal/System/Time.h"
#include "klee/Internal/Support/ErrorHandling.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/CommandLine.h"

#if LLVM_VERSION_CODE < LLVM_VERSION(3, 5)
#include "llvm/Support/CallSite.h"
#else
#include "llvm/IR/CallSite.h"
#endif

#include <cassert>
#include <fstream>
#include <climits>

using namespace klee;
using namespace llvm;

namespace {
  cl::opt<bool>
  DebugLogMerge("debug-log-merge");
}

namespace klee {
  extern RNG theRNG;
}

Searcher::~Searcher() {
}

/** 
	#####################
  # ValueSearch Start #
	#####################
*/
#include "llvm/Support/Errno.h"
#include "klee/Constraints.h"
#include "klee/util/ExprPPrinter.h"
#include <map>
#include <string>
class FeatureExtractor {
	public: 
	
		/**
		Summarize ExecutionState into a feature vector
		*/
		static void extract(ExecutionState * current, std::map<std::string, double> &dbl_map, std::map<std::string, std::string> &str_map, int debug) {
			
			// Instructions Features
			if (true) { 
				// get type name
				std::string pc_type_str;
				llvm::raw_string_ostream pc_rso(pc_type_str);
				current->pc->inst->getType()->print(pc_rso);
				dbl_map.insert(std::pair<std::string, double>("instr_pc_getNumOperands", (double) current->pc->inst->getNumOperands()));
				dbl_map.insert(std::pair<std::string, double>("instr_pc_getOpcode", (double) current->pc->inst->getOpcode()));
				dbl_map.insert(std::pair<std::string, double>("instr_pc_use_empty", (double) current->pc->inst->use_empty()));
				str_map.insert(std::pair<std::string, std::string>("instr_pc_getType", pc_rso.str().c_str()));
				
				std::string prevpc_type_str;
				llvm::raw_string_ostream prevpc_rso(prevpc_type_str);
				current->prevPC->inst->getType()->print(prevpc_rso);
				dbl_map.insert(std::pair<std::string, double>("instr_prevpc_getNumOperands", (double) current->prevPC->inst->getNumOperands()));
				dbl_map.insert(std::pair<std::string, double>("instr_prevpc_getOpcode", (double) current->prevPC->inst->getOpcode()));
				dbl_map.insert(std::pair<std::string, double>("instr_prevpc_use_empty", (double) current->prevPC->inst->use_empty()));
				str_map.insert(std::pair<std::string, std::string>("instr_prevpc_getType", prevpc_rso.str().c_str()));
			}
			
			// WeightedRandomSearcher features
			if (true) {		
				dbl_map.insert(std::pair<std::string, double>("weights_Depth", get_random_searcher_weights(current, "Depth")));
				dbl_map.insert(std::pair<std::string, double>("weights_InstCount", get_random_searcher_weights(current, "InstCount")));
				dbl_map.insert(std::pair<std::string, double>("weights_CPInstCount", get_random_searcher_weights(current, "CPInstCount")));
				dbl_map.insert(std::pair<std::string, double>("weights_QueryCost", get_random_searcher_weights(current, "QueryCost")));
				dbl_map.insert(std::pair<std::string, double>("weights_CoveringNew", get_random_searcher_weights(current, "CoveringNew")));
				dbl_map.insert(std::pair<std::string, double>("weights_MinDistToUncovered", get_random_searcher_weights(current, "MinDistToUncovered")));
				dbl_map.insert(std::pair<std::string, double>("weights_coveredNew", (double) current->coveredNew));
				dbl_map.insert(std::pair<std::string, double>("weights_forkDisabled", (double) current->forkDisabled));
			}
			
			// Constraints Features
			if (true) {
				
				// get constraints string representation
				const char * constr = get_constraints(current);
				
				// Propositional SAT features
				
					// suggested by "Hardness Estimation of QFBV SMT Problems"	
					int num_ands = count_occurences(constr, (char *)"And");
					int num_asserted = (int) current->constraints.size(); //# of asserted formulas
					int num_eq = count_occurences(constr, (char *)"Eq");
					int num_xor = count_occurences(constr, (char *)"Xor");
					int num_bool = count_occurences(constr, (char *)"w"); //num of boolean varialbes to be assigned
					int num_theory = count_occurences(constr, (char *)"("); //# of theory athoms
					int clause_to_vars = (num_asserted + num_ands + 2 * (num_eq - 2) + num_xor) / (num_bool + num_theory);
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_ands", (double) num_ands ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_asserted", (double) num_asserted ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_eq", (double) num_eq ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_xor", (double) num_xor ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_bool", (double) num_bool ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_num_theory", (double) num_theory ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_clause_to_vars", (double) clause_to_vars ));
					dbl_map.insert(std::pair<std::string, double>("sat_paper_vars_clause", (double) (clause_to_vars == 0) ? 1 : (1 / clause_to_vars) ));

					// original
				
					// arithmetic
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_sub", (double) count_occurences(constr, (char *)"Sub") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_mul", (double) count_occurences(constr, (char *)"Mul") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_udiv", (double) count_occurences(constr, (char *)"UDiv") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_sdiv", (double) count_occurences(constr, (char *)"SDiv") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_urem", (double) count_occurences(constr, (char *)"URem") ));
					
					// bitwise
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_or", (double) count_occurences(constr, (char *)"Or") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_shl", (double) count_occurences(constr, (char *)"Shl") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_lshr", (double) count_occurences(constr, (char *)"LShr") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_ashr", (double) count_occurences(constr, (char *)"AShr") ));

					// Bitvector Manipulation
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_extract", (double) count_occurences(constr, (char *)"Extract") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_concat", (double) count_occurences(constr, (char *)"Concat") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_zext", (double) count_occurences(constr, (char *)"ZExt") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_sext", (double) count_occurences(constr, (char *)"SExt") ));
					
					// Macros
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_lsb", (double) count_occurences(constr, (char *)"ReadLSB") ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_msb", (double) count_occurences(constr, (char *)"ReadMSB") ));
					
				// Propositional SAT features
					double clause_depth [3];
					consecutive_occurences_stats(constr, '(', ')', clause_depth); // {std, avg, max}
					dbl_map.insert(std::pair<std::string, double>("sat_org_cls_std", clause_depth[0] ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_cls_avg", clause_depth[1] ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_cls_max", clause_depth[2] ));
					dbl_map.insert(std::pair<std::string, double>("sat_org_num_capitals", (double) count_capitals(constr) ));
			}
			
			
			// Debug
			if (debug) {
				klee_message("####################################");
				klee_message("######### State Features ###########");
				for(std::map<std::string, double>::iterator iter = dbl_map.begin(); iter != dbl_map.end(); ++iter) {
					klee_message("%s=%.2f", iter->first.c_str(), iter->second);
				}
				for(std::map<std::string, std::string>::iterator iter = str_map.begin(); iter != str_map.end(); ++iter) {
					klee_message("%s=%s", iter->first.c_str(), iter->second.c_str());
				}
				klee_message("####################################");
			}
		}
		
	private:

		/* String-related methods */
		
		/**
		Count occurences of substr in str
		*/
		static int count_occurences(const char * str, char * substr) {
			int count = 0;
			while((str = strstr(str, substr))) {
				 count++;
				 str++;
			}
			return count;
		}
		
		/**
		Count the number of capital letters in a string
		*/
		static int count_capitals(const char * str) {
				int count = 0;
				int i = 0;
				while (str[i] != '\0') {
					if (isupper(str[i])) {count++;}
					i++;
				}
				return count;   
		}
		
		/**
		Given a String str, a search_char and a break char, calculates the 
		length of all possible sequences of consequtive search_char's without a middle break_char.
		These possible sequences of length are aggregated into the std, avg and max length, 
		and saved into the first three places in results
		
		@assert len(results) >= 3
		*/
		static void consecutive_occurences_stats(const char * str, char search_char, char break_char, double * results) {
			results[0] = 1.0; // std
			results[1] = 0.0; // avg
			results[2] = 0.0; // max
			
			// helpers
			int cur_mean, delta1, delta2;
			
			int current = 0;
			int count_switches = 0;
			int ind = 0;
			while(true) {
					if (str[ind] == search_char) {
						current++;
					} else if (str[ind] == break_char or str[ind] == '\0') {
						
						// update n
						count_switches++;
						
						// update max
						results[2] = (results[2] < (double) current) ? (double) current : results[2];
						
						// online update std
						cur_mean = (results[1] / (count_switches - 1));
						delta1 = current - cur_mean;
						delta2 = current - cur_mean - (delta1 / count_switches);
						results[0] += delta1*delta2;  // update sum of deltas^2
						
						// update sum (which will become avg)
						results[1] = results[1] + (double) current;
					
						// update current for next one
						current--;
						if (str[ind] == '\0') { break; }
				 }
				 ind++;
			}
			results[0] = results[0] / (count_switches - 1);  // sum of deltas^2 -> std
			results[1] = (count_switches == 0) ? 0.0 : results[1] / (double) count_switches; // sum -> avg
		}
		
		
		/* ExecutionState-related methods */
		
		/**
		Extracts a string representation in KQuery format of the current constraints, from the Execution state
		*/
		static const char * get_constraints(ExecutionState * state) {
			std::string Str;
			llvm::raw_string_ostream info(Str);
			ExprPPrinter::printConstraints(info, state->constraints);
			std::string res = info.str();

			return res.c_str();
		}
		
		/**
		Given an ExecutionState and a choosen WeightType in string format,
		Extract the weight given by WeightedRandomSearcher to this state with this weight type.
		Possible types are: Depth, InstCount, CPInstCount, QueryCost, CoveringNew, MinDistToUncovered
		*/
		static double get_random_searcher_weights(ExecutionState *es, const std::string type) {
			if (type == "Depth") {
					return es->weight;
			} else if (type == "InstCount") {
					uint64_t count = theStatisticManager->getIndexedValue(stats::instructions,
																															es->pc->info->id);
					double inv = 1. / std::max((uint64_t) 1, count);
					return inv * inv;
			} else if (type == "CPInstCount") {
					StackFrame &sf = es->stack.back();
					uint64_t count = sf.callPathNode->statistics.getValue(stats::instructions);
					double inv = 1. / std::max((uint64_t) 1, count);
					return inv;
			} else if (type == "QueryCost") {
					return (es->queryCost < .1) ? 1. : 1./es->queryCost;
			} else if (type == "CoveringNew") {
					uint64_t md2u = computeMinDistToUncovered(es->pc,
																									es->stack.back().minDistToUncoveredOnReturn);
					double invMD2U = 1. / (md2u ? md2u : 10000);
					double invCovNew = 0.;
					if (es->instsSinceCovNew)
						invCovNew = 1. / std::max(1, (int) es->instsSinceCovNew - 1000);
					return (invCovNew * invCovNew + invMD2U * invMD2U);
			} else if (type == "MinDistToUncovered") {
					uint64_t md2u = computeMinDistToUncovered(es->pc,
																									es->stack.back().minDistToUncoveredOnReturn);
					double invMD2U = 1. / (md2u ? md2u : 10000);
					return invMD2U * invMD2U;
			}
			return -1.0;
		}
};

/** 
	###################
  # ValueSearch End #
	###################
*/

///

ExecutionState &DFSSearcher::selectState() {
  return *states.back();
}

void DFSSearcher::update(ExecutionState *current,
                         const std::vector<ExecutionState *> &addedStates,
                         const std::vector<ExecutionState *> &removedStates) {
  states.insert(states.end(),
                addedStates.begin(),
                addedStates.end());
  for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
                                                     ie = removedStates.end();
       it != ie; ++it) {
    ExecutionState *es = *it;
    if (es == states.back()) {
      states.pop_back();
    } else {
      bool ok = false;

      for (std::vector<ExecutionState*>::iterator it = states.begin(),
             ie = states.end(); it != ie; ++it) {
        if (es==*it) {
          states.erase(it);
          ok = true;
          break;
        }
      }

      (void) ok;
      assert(ok && "invalid state removed");
    }
  }
}

///

ExecutionState &BFSSearcher::selectState() {
	/** ValueSearch - Usage Example Start */
	std::map<std::string, double> dbl_map;
	std::map<std::string, std::string> str_map;
	FeatureExtractor::extract(&(*states.front()), dbl_map, str_map, 1);
	/** ValueSearch - Usage Example End */  
  return *states.front();
}

void BFSSearcher::update(ExecutionState *current,
                         const std::vector<ExecutionState *> &addedStates,
                         const std::vector<ExecutionState *> &removedStates) {  
  // Assumption: If new states were added KLEE forked, therefore states evolved.
  // constraints were added to the current state, it evolved.
  if (!addedStates.empty() && current &&
      std::find(removedStates.begin(), removedStates.end(), current) ==
          removedStates.end()) {
    auto pos = std::find(states.begin(), states.end(), current);
    assert(pos != states.end());
    states.erase(pos);
    states.push_back(current);
  }

  states.insert(states.end(),
                addedStates.begin(),
                addedStates.end());
  for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
                                                     ie = removedStates.end();
       it != ie; ++it) {
    ExecutionState *es = *it;
    if (es == states.front()) {
      states.pop_front();
    } else {
      bool ok = false;

      for (std::deque<ExecutionState*>::iterator it = states.begin(),
             ie = states.end(); it != ie; ++it) {
        if (es==*it) {
          states.erase(it);
          ok = true;
          break;
        }
      }

      (void) ok;
      assert(ok && "invalid state removed");
    }
  }
}

///

ExecutionState &RandomSearcher::selectState() {
  return *states[theRNG.getInt32()%states.size()];
}

void
RandomSearcher::update(ExecutionState *current,
                       const std::vector<ExecutionState *> &addedStates,
                       const std::vector<ExecutionState *> &removedStates) {
  states.insert(states.end(),
                addedStates.begin(),
                addedStates.end());
  for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
                                                     ie = removedStates.end();
       it != ie; ++it) {
    ExecutionState *es = *it;
    bool ok = false;

    for (std::vector<ExecutionState*>::iterator it = states.begin(),
           ie = states.end(); it != ie; ++it) {
      if (es==*it) {
        states.erase(it);
        ok = true;
        break;
      }
    }
    
    assert(ok && "invalid state removed");
  }
}

///

WeightedRandomSearcher::WeightedRandomSearcher(WeightType _type)
  : states(new DiscretePDF<ExecutionState*>()),
    type(_type) {
  switch(type) {
  case Depth: 
    updateWeights = false;
    break;
  case InstCount:
  case CPInstCount:
  case QueryCost:
  case MinDistToUncovered:
  case CoveringNew:
    updateWeights = true;
    break;
  default:
    assert(0 && "invalid weight type");
  }
}

WeightedRandomSearcher::~WeightedRandomSearcher() {
  delete states;
}

ExecutionState &WeightedRandomSearcher::selectState() {
  return *states->choose(theRNG.getDoubleL());
}

double WeightedRandomSearcher::getWeight(ExecutionState *es) {
  switch(type) {
  default:
  case Depth: 
    return es->weight;
  case InstCount: {
    uint64_t count = theStatisticManager->getIndexedValue(stats::instructions,
                                                          es->pc->info->id);
    double inv = 1. / std::max((uint64_t) 1, count);
    return inv * inv;
  }
  case CPInstCount: {
    StackFrame &sf = es->stack.back();
    uint64_t count = sf.callPathNode->statistics.getValue(stats::instructions);
    double inv = 1. / std::max((uint64_t) 1, count);
    return inv;
  }
  case QueryCost:
    return (es->queryCost < .1) ? 1. : 1./es->queryCost;
  case CoveringNew:
  case MinDistToUncovered: {
    uint64_t md2u = computeMinDistToUncovered(es->pc,
                                              es->stack.back().minDistToUncoveredOnReturn);

    double invMD2U = 1. / (md2u ? md2u : 10000);
    if (type==CoveringNew) {
      double invCovNew = 0.;
      if (es->instsSinceCovNew)
        invCovNew = 1. / std::max(1, (int) es->instsSinceCovNew - 1000);
      return (invCovNew * invCovNew + invMD2U * invMD2U);
    } else {
      return invMD2U * invMD2U;
    }
  }
  }
}

void WeightedRandomSearcher::update(
    ExecutionState *current, const std::vector<ExecutionState *> &addedStates,
    const std::vector<ExecutionState *> &removedStates) {
  if (current && updateWeights &&
      std::find(removedStates.begin(), removedStates.end(), current) ==
          removedStates.end())
    states->update(current, getWeight(current));

  for (std::vector<ExecutionState *>::const_iterator it = addedStates.begin(),
                                                     ie = addedStates.end();
       it != ie; ++it) {
    ExecutionState *es = *it;
    states->insert(es, getWeight(es));
  }

  for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
                                                     ie = removedStates.end();
       it != ie; ++it) {
    states->remove(*it);
  }
}

bool WeightedRandomSearcher::empty() { 
  return states->empty(); 
}

///

RandomPathSearcher::RandomPathSearcher(Executor &_executor)
  : executor(_executor) {
}

RandomPathSearcher::~RandomPathSearcher() {
}

ExecutionState &RandomPathSearcher::selectState() {
  unsigned flips=0, bits=0;
  PTree::Node *n = executor.processTree->root;
  
  while (!n->data) {
    if (!n->left) {
      n = n->right;
    } else if (!n->right) {
      n = n->left;
    } else {
      if (bits==0) {
        flips = theRNG.getInt32();
        bits = 32;
      }
      --bits;
      n = (flips&(1<<bits)) ? n->left : n->right;
    }
  }

  return *n->data;
}

void
RandomPathSearcher::update(ExecutionState *current,
                           const std::vector<ExecutionState *> &addedStates,
                           const std::vector<ExecutionState *> &removedStates) {
}

bool RandomPathSearcher::empty() { 
  return executor.states.empty(); 
}

///

BatchingSearcher::BatchingSearcher(Searcher *_baseSearcher,
                                   double _timeBudget,
                                   unsigned _instructionBudget) 
  : baseSearcher(_baseSearcher),
    timeBudget(_timeBudget),
    instructionBudget(_instructionBudget),
    lastState(0) {
  
}

BatchingSearcher::~BatchingSearcher() {
  delete baseSearcher;
}

ExecutionState &BatchingSearcher::selectState() {
  if (!lastState || 
      (util::getWallTime()-lastStartTime)>timeBudget ||
      (stats::instructions-lastStartInstructions)>instructionBudget) {
    if (lastState) {
      double delta = util::getWallTime()-lastStartTime;
      if (delta>timeBudget*1.1) {
        klee_message("increased time budget from %f to %f\n", timeBudget,
                     delta);
        timeBudget = delta;
      }
    }
    lastState = &baseSearcher->selectState();
    lastStartTime = util::getWallTime();
    lastStartInstructions = stats::instructions;
    return *lastState;
  } else {
    return *lastState;
  }
}

void
BatchingSearcher::update(ExecutionState *current,
                         const std::vector<ExecutionState *> &addedStates,
                         const std::vector<ExecutionState *> &removedStates) {
  if (std::find(removedStates.begin(), removedStates.end(), lastState) !=
      removedStates.end())
    lastState = 0;
  baseSearcher->update(current, addedStates, removedStates);
}

/***/

IterativeDeepeningTimeSearcher::IterativeDeepeningTimeSearcher(Searcher *_baseSearcher)
  : baseSearcher(_baseSearcher),
    time(1.) {
}

IterativeDeepeningTimeSearcher::~IterativeDeepeningTimeSearcher() {
  delete baseSearcher;
}

ExecutionState &IterativeDeepeningTimeSearcher::selectState() {
  ExecutionState &res = baseSearcher->selectState();
  startTime = util::getWallTime();
  return res;
}

void IterativeDeepeningTimeSearcher::update(
    ExecutionState *current, const std::vector<ExecutionState *> &addedStates,
    const std::vector<ExecutionState *> &removedStates) {
  double elapsed = util::getWallTime() - startTime;

  if (!removedStates.empty()) {
    std::vector<ExecutionState *> alt = removedStates;
    for (std::vector<ExecutionState *>::const_iterator
             it = removedStates.begin(),
             ie = removedStates.end();
         it != ie; ++it) {
      ExecutionState *es = *it;
      std::set<ExecutionState*>::const_iterator it2 = pausedStates.find(es);
      if (it2 != pausedStates.end()) {
        pausedStates.erase(it2);
        alt.erase(std::remove(alt.begin(), alt.end(), es), alt.end());
      }
    }    
    baseSearcher->update(current, addedStates, alt);
  } else {
    baseSearcher->update(current, addedStates, removedStates);
  }

  if (current &&
      std::find(removedStates.begin(), removedStates.end(), current) ==
          removedStates.end() &&
      elapsed > time) {
    pausedStates.insert(current);
    baseSearcher->removeState(current);
  }

  if (baseSearcher->empty()) {
    time *= 2;
    klee_message("increased time budget to %f\n", time);
    std::vector<ExecutionState *> ps(pausedStates.begin(), pausedStates.end());
    baseSearcher->update(0, ps, std::vector<ExecutionState *>());
    pausedStates.clear();
  }
}

/***/

InterleavedSearcher::InterleavedSearcher(const std::vector<Searcher*> &_searchers)
  : searchers(_searchers),
    index(1) {
}

InterleavedSearcher::~InterleavedSearcher() {
  for (std::vector<Searcher*>::const_iterator it = searchers.begin(),
         ie = searchers.end(); it != ie; ++it)
    delete *it;
}

ExecutionState &InterleavedSearcher::selectState() {
  Searcher *s = searchers[--index];
  if (index==0) index = searchers.size();
  return s->selectState();
}

void InterleavedSearcher::update(
    ExecutionState *current, const std::vector<ExecutionState *> &addedStates,
    const std::vector<ExecutionState *> &removedStates) {
  for (std::vector<Searcher*>::const_iterator it = searchers.begin(),
         ie = searchers.end(); it != ie; ++it)
    (*it)->update(current, addedStates, removedStates);
}
