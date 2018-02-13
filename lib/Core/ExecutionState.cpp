//===-- ExecutionState.cpp ------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "klee/ExecutionState.h"

#include "klee/Internal/Module/Cell.h"
#include "klee/Internal/Module/InstructionInfoTable.h"
#include "klee/Internal/Module/KInstruction.h"
#include "klee/Internal/Module/KModule.h"
#include "klee/Expr.h"

#include "Memory.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/raw_ostream.h"

#include <iomanip>
#include <sstream>
#include <cassert>
#include <map>
#include <set>
#include <stdarg.h>

//@COMPILATION_PROJECT:
#include "klee/util/StringUtils.h"
#include "klee/util/ExprPPrinter.h"
#include "klee/Statistics.h"
#include "llvm/Support/Errno.h"
#include "klee/Constraints.h"
#include "StatsTracker.h"
#include "CoreStats.h"
#include "klee/Internal/Support/ErrorHandling.h"

using namespace llvm;
using namespace klee;

namespace { 
  cl::opt<bool>
  DebugLogStateMerge("debug-log-state-merge");
}

/***/

StackFrame::StackFrame(KInstIterator _caller, KFunction *_kf)
  : caller(_caller), kf(_kf), callPathNode(0), 
    minDistToUncoveredOnReturn(0), varargs(0) {
  locals = new Cell[kf->numRegisters];
}

StackFrame::StackFrame(const StackFrame &s) 
  : caller(s.caller),
    kf(s.kf),
    callPathNode(s.callPathNode),
    allocas(s.allocas),
    minDistToUncoveredOnReturn(s.minDistToUncoveredOnReturn),
    varargs(s.varargs) {
  locals = new Cell[s.kf->numRegisters];
  for (unsigned i=0; i<s.kf->numRegisters; i++)
    locals[i] = s.locals[i];
}

StackFrame::~StackFrame() { 
  delete[] locals; 
}

/***/

ExecutionState::ExecutionState(KFunction *kf) :
    pc(kf->instructions),
    prevPC(pc),

    queryCost(0.), 
    weight(1),
    depth(0),

    instsSinceCovNew(0),
    coveredNew(false),
    forkDisabled(false),
    ptreeNode(0) {
  pushFrame(0, kf);
}

ExecutionState::ExecutionState(const std::vector<ref<Expr> > &assumptions)
    : constraints(assumptions), queryCost(0.), ptreeNode(0) {}

ExecutionState::~ExecutionState() {
  for (unsigned int i=0; i<symbolics.size(); i++)
  {
    const MemoryObject *mo = symbolics[i].first;
    assert(mo->refCount > 0);
    mo->refCount--;
    if (mo->refCount == 0)
      delete mo;
  }

  while (!stack.empty()) popFrame();
}

ExecutionState::ExecutionState(const ExecutionState& state):
    fnAliases(state.fnAliases),
    pc(state.pc),
    prevPC(state.prevPC),
    stack(state.stack),
    incomingBBIndex(state.incomingBBIndex),

    addressSpace(state.addressSpace),
    constraints(state.constraints),

    queryCost(state.queryCost),
    weight(state.weight),
    depth(state.depth),

    pathOS(state.pathOS),
    symPathOS(state.symPathOS),

    instsSinceCovNew(state.instsSinceCovNew),
    coveredNew(state.coveredNew),
    forkDisabled(state.forkDisabled),
    coveredLines(state.coveredLines),
    ptreeNode(state.ptreeNode),
    symbolics(state.symbolics),
    arrayNames(state.arrayNames)
{
  for (unsigned int i=0; i<symbolics.size(); i++)
    symbolics[i].first->refCount++;
}

ExecutionState *ExecutionState::branch() {
  depth++;

  ExecutionState *falseState = new ExecutionState(*this);
  falseState->coveredNew = false;
  falseState->coveredLines.clear();

  weight *= .5;
  falseState->weight -= weight;

  return falseState;
}

void ExecutionState::pushFrame(KInstIterator caller, KFunction *kf) {
  stack.push_back(StackFrame(caller,kf));
}

void ExecutionState::popFrame() {
  StackFrame &sf = stack.back();
  for (std::vector<const MemoryObject*>::iterator it = sf.allocas.begin(), 
         ie = sf.allocas.end(); it != ie; ++it)
    addressSpace.unbindObject(*it);
  stack.pop_back();
}

void ExecutionState::addSymbolic(const MemoryObject *mo, const Array *array) { 
  mo->refCount++;
  symbolics.push_back(std::make_pair(mo, array));
}
///

std::string ExecutionState::getFnAlias(std::string fn) {
  std::map < std::string, std::string >::iterator it = fnAliases.find(fn);
  if (it != fnAliases.end())
    return it->second;
  else return "";
}

void ExecutionState::addFnAlias(std::string old_fn, std::string new_fn) {
  fnAliases[old_fn] = new_fn;
}

void ExecutionState::removeFnAlias(std::string fn) {
  fnAliases.erase(fn);
}

/**/

llvm::raw_ostream &klee::operator<<(llvm::raw_ostream &os, const MemoryMap &mm) {
  os << "{";
  MemoryMap::iterator it = mm.begin();
  MemoryMap::iterator ie = mm.end();
  if (it!=ie) {
    os << "MO" << it->first->id << ":" << it->second;
    for (++it; it!=ie; ++it)
      os << ", MO" << it->first->id << ":" << it->second;
  }
  os << "}";
  return os;
}

bool ExecutionState::merge(const ExecutionState &b) {
  if (DebugLogStateMerge)
    llvm::errs() << "-- attempting merge of A:" << this << " with B:" << &b
                 << "--\n";
  if (pc != b.pc)
    return false;

  // XXX is it even possible for these to differ? does it matter? probably
  // implies difference in object states?
  if (symbolics!=b.symbolics)
    return false;

  {
    std::vector<StackFrame>::const_iterator itA = stack.begin();
    std::vector<StackFrame>::const_iterator itB = b.stack.begin();
    while (itA!=stack.end() && itB!=b.stack.end()) {
      // XXX vaargs?
      if (itA->caller!=itB->caller || itA->kf!=itB->kf)
        return false;
      ++itA;
      ++itB;
    }
    if (itA!=stack.end() || itB!=b.stack.end())
      return false;
  }

  std::set< ref<Expr> > aConstraints(constraints.begin(), constraints.end());
  std::set< ref<Expr> > bConstraints(b.constraints.begin(), 
                                     b.constraints.end());
  std::set< ref<Expr> > commonConstraints, aSuffix, bSuffix;
  std::set_intersection(aConstraints.begin(), aConstraints.end(),
                        bConstraints.begin(), bConstraints.end(),
                        std::inserter(commonConstraints, commonConstraints.begin()));
  std::set_difference(aConstraints.begin(), aConstraints.end(),
                      commonConstraints.begin(), commonConstraints.end(),
                      std::inserter(aSuffix, aSuffix.end()));
  std::set_difference(bConstraints.begin(), bConstraints.end(),
                      commonConstraints.begin(), commonConstraints.end(),
                      std::inserter(bSuffix, bSuffix.end()));
  if (DebugLogStateMerge) {
    llvm::errs() << "\tconstraint prefix: [";
    for (std::set<ref<Expr> >::iterator it = commonConstraints.begin(),
                                        ie = commonConstraints.end();
         it != ie; ++it)
      llvm::errs() << *it << ", ";
    llvm::errs() << "]\n";
    llvm::errs() << "\tA suffix: [";
    for (std::set<ref<Expr> >::iterator it = aSuffix.begin(),
                                        ie = aSuffix.end();
         it != ie; ++it)
      llvm::errs() << *it << ", ";
    llvm::errs() << "]\n";
    llvm::errs() << "\tB suffix: [";
    for (std::set<ref<Expr> >::iterator it = bSuffix.begin(),
                                        ie = bSuffix.end();
         it != ie; ++it)
      llvm::errs() << *it << ", ";
    llvm::errs() << "]\n";
  }

  // We cannot merge if addresses would resolve differently in the
  // states. This means:
  // 
  // 1. Any objects created since the branch in either object must
  // have been free'd.
  //
  // 2. We cannot have free'd any pre-existing object in one state
  // and not the other

  if (DebugLogStateMerge) {
    llvm::errs() << "\tchecking object states\n";
    llvm::errs() << "A: " << addressSpace.objects << "\n";
    llvm::errs() << "B: " << b.addressSpace.objects << "\n";
  }
    
  std::set<const MemoryObject*> mutated;
  MemoryMap::iterator ai = addressSpace.objects.begin();
  MemoryMap::iterator bi = b.addressSpace.objects.begin();
  MemoryMap::iterator ae = addressSpace.objects.end();
  MemoryMap::iterator be = b.addressSpace.objects.end();
  for (; ai!=ae && bi!=be; ++ai, ++bi) {
    if (ai->first != bi->first) {
      if (DebugLogStateMerge) {
        if (ai->first < bi->first) {
          llvm::errs() << "\t\tB misses binding for: " << ai->first->id << "\n";
        } else {
          llvm::errs() << "\t\tA misses binding for: " << bi->first->id << "\n";
        }
      }
      return false;
    }
    if (ai->second != bi->second) {
      if (DebugLogStateMerge)
        llvm::errs() << "\t\tmutated: " << ai->first->id << "\n";
      mutated.insert(ai->first);
    }
  }
  if (ai!=ae || bi!=be) {
    if (DebugLogStateMerge)
      llvm::errs() << "\t\tmappings differ\n";
    return false;
  }
  
  // merge stack

  ref<Expr> inA = ConstantExpr::alloc(1, Expr::Bool);
  ref<Expr> inB = ConstantExpr::alloc(1, Expr::Bool);
  for (std::set< ref<Expr> >::iterator it = aSuffix.begin(), 
         ie = aSuffix.end(); it != ie; ++it)
    inA = AndExpr::create(inA, *it);
  for (std::set< ref<Expr> >::iterator it = bSuffix.begin(), 
         ie = bSuffix.end(); it != ie; ++it)
    inB = AndExpr::create(inB, *it);

  // XXX should we have a preference as to which predicate to use?
  // it seems like it can make a difference, even though logically
  // they must contradict each other and so inA => !inB

  std::vector<StackFrame>::iterator itA = stack.begin();
  std::vector<StackFrame>::const_iterator itB = b.stack.begin();
  for (; itA!=stack.end(); ++itA, ++itB) {
    StackFrame &af = *itA;
    const StackFrame &bf = *itB;
    for (unsigned i=0; i<af.kf->numRegisters; i++) {
      ref<Expr> &av = af.locals[i].value;
      const ref<Expr> &bv = bf.locals[i].value;
      if (av.isNull() || bv.isNull()) {
        // if one is null then by implication (we are at same pc)
        // we cannot reuse this local, so just ignore
      } else {
        av = SelectExpr::create(inA, av, bv);
      }
    }
  }

  for (std::set<const MemoryObject*>::iterator it = mutated.begin(), 
         ie = mutated.end(); it != ie; ++it) {
    const MemoryObject *mo = *it;
    const ObjectState *os = addressSpace.findObject(mo);
    const ObjectState *otherOS = b.addressSpace.findObject(mo);
    assert(os && !os->readOnly && 
           "objects mutated but not writable in merging state");
    assert(otherOS);

    ObjectState *wos = addressSpace.getWriteable(mo, os);
    for (unsigned i=0; i<mo->size; i++) {
      ref<Expr> av = wos->read8(i);
      ref<Expr> bv = otherOS->read8(i);
      wos->write(i, SelectExpr::create(inA, av, bv));
    }
  }

  constraints = ConstraintManager();
  for (std::set< ref<Expr> >::iterator it = commonConstraints.begin(), 
         ie = commonConstraints.end(); it != ie; ++it)
    constraints.addConstraint(*it);
  constraints.addConstraint(OrExpr::create(inA, inB));

  return true;
}

void ExecutionState::dumpStack(llvm::raw_ostream &out) const {
  unsigned idx = 0;
  const KInstruction *target = prevPC;
  for (ExecutionState::stack_ty::const_reverse_iterator
         it = stack.rbegin(), ie = stack.rend();
       it != ie; ++it) {
    const StackFrame &sf = *it;
    Function *f = sf.kf->function;
    const InstructionInfo &ii = *target->info;
    out << "\t#" << idx++;
    std::stringstream AssStream;
    AssStream << std::setw(8) << std::setfill('0') << ii.assemblyLine;
    out << AssStream.str();
    out << " in " << f->getName().str() << " (";
    // Yawn, we could go up and print varargs if we wanted to.
    unsigned index = 0;
    for (Function::arg_iterator ai = f->arg_begin(), ae = f->arg_end();
         ai != ae; ++ai) {
      if (ai!=f->arg_begin()) out << ", ";

      out << ai->getName().str();
      // XXX should go through function
      ref<Expr> value = sf.locals[sf.kf->getArgRegister(index++)].value;
      if (value.get() && isa<ConstantExpr>(value))
        out << "=" << value;
    }
    out << ")";
    if (ii.file != "")
      out << " at " << ii.file << ":" << ii.line;
    out << "\n";
    target = sf.caller;
  }
}


//@COMPILATION_PROJECT
void ExecutionState::extractStateMaps(std::map<std::string, double> &dbl_map, std::map<std::string, std::string> &str_map, bool debug){
  // Instructions Features
    // get type name
    std::string pc_type_str;
    llvm::raw_string_ostream pc_rso(pc_type_str);
    this->pc->inst->getType()->print(pc_rso);
    dbl_map.insert(std::pair<std::string, double>(INSTR_PC_OPERANDS, (double) this->pc->inst->getNumOperands()));
    dbl_map.insert(std::pair<std::string, double>(INSTR_PC_OPCODE, (double) this->pc->inst->getOpcode()));
    dbl_map.insert(std::pair<std::string, double>(INSTR_PC_USE_EMPTY, (double) this->pc->inst->use_empty()));
    str_map.insert(std::pair<std::string, std::string>(INSTR_PC_TYPE, pc_rso.str().c_str()));

    std::string prevpc_type_str;
    llvm::raw_string_ostream prevpc_rso(prevpc_type_str);
    this->prevPC->inst->getType()->print(prevpc_rso);
    dbl_map.insert(std::pair<std::string, double>(INSTR_PPC_OPERANDS, (double) this->prevPC->inst->getNumOperands()));
    dbl_map.insert(std::pair<std::string, double>(INSTR_PPC_OPCODE, (double) this->prevPC->inst->getOpcode()));
    dbl_map.insert(std::pair<std::string, double>(INSTR_PPC_USE_EMPTY, (double) this->prevPC->inst->use_empty()));
    str_map.insert(std::pair<std::string, std::string>(INSTR_PPC_TYPE, prevpc_rso.str().c_str()));

    // Constraints Features
    // get constraints string representation
    const char * constr = get_constraints();

    // Propositional SAT features

    // suggested by "Hardness Estimation of QFBV SMT Problems"
    int num_ands = StringUtils::count_occurences(constr, (char *)"And");
    int num_asserted = (int) this->constraints.size(); //# of asserted formulas
    int num_eq = StringUtils::count_occurences(constr, (char *)"Eq");
    int num_xor = StringUtils::count_occurences(constr, (char *)"Xor");
    int num_bool = StringUtils::count_occurences(constr, (char *)"w"); //num of boolean varialbes to be assigned
    int num_theory = StringUtils::count_occurences(constr, (char *)"("); //# of theory athoms
    int clause_to_vars = (num_asserted + num_ands + 2 * (num_eq - 2) + num_xor) / (num_bool + num_theory);
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_ANDS, (double) num_ands ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_ASSERTED, (double) num_asserted ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_EQUAL, (double) num_eq ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_XOR, (double) num_xor ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_BOOL, (double) num_bool ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_THEORY, (double) num_theory ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_CLAUSES_VARS, (double) clause_to_vars ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_VARS_CLAUSES, (double) (clause_to_vars == 0) ? 1 : (1 / clause_to_vars) ));
    // original
    // arithmetic
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_SUB, (double) StringUtils::count_occurences(constr, (char *)"Sub") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_MUL, (double) StringUtils::count_occurences(constr, (char *)"Mul") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_UDIV, (double) StringUtils::count_occurences(constr, (char *)"UDiv") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_SDIV, (double) StringUtils::count_occurences(constr, (char *)"SDiv") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_UREM, (double) StringUtils::count_occurences(constr, (char *)"URem") ));
    // bitwise
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_OR, (double) StringUtils::count_occurences(constr, (char *)"Or") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_SHL, (double) StringUtils::count_occurences(constr, (char *)"Shl") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_LSHR, (double) StringUtils::count_occurences(constr, (char *)"LShr") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_ASHR, (double) StringUtils::count_occurences(constr, (char *)"AShr") ));
    // Bitvector Manipulation
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_EXTRACT, (double) StringUtils::count_occurences(constr, (char *)"Extract") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_CONCAT, (double) StringUtils::count_occurences(constr, (char *)"Concat") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_ZEXT, (double) StringUtils::count_occurences(constr, (char *)"ZExt") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_SEXT, (double) StringUtils::count_occurences(constr, (char *)"SExt") ));

    // Macros
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_LSB, (double) StringUtils::count_occurences(constr, (char *)"ReadLSB") ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_MSB, (double) StringUtils::count_occurences(constr, (char *)"ReadMSB") ));

    // Propositional SAT features
    double clause_depth [3];
    StringUtils::consecutive_occurences_stats(constr, '(', ')', clause_depth); // {std, avg, max}
    dbl_map.insert(std::pair<std::string, double>(SAT_CLS_STD, clause_depth[0] ));
    dbl_map.insert(std::pair<std::string, double>(SAT_CLS_AVG, clause_depth[1] ));
    dbl_map.insert(std::pair<std::string, double>(SAT_CLS_MAX, clause_depth[2] ));
    dbl_map.insert(std::pair<std::string, double>(SAT_NUM_CAPS, (double) StringUtils::count_capitals(constr) ));
    // WeightedRandomSearcher features
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_DEPTH, get_random_searcher_weights("Depth")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_INST_CNT, get_random_searcher_weights("InstCount")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_CPI_CNT, get_random_searcher_weights("CPInstCount")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_QUERY_CST, get_random_searcher_weights("QueryCost")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_COVER_NEW, get_random_searcher_weights("CoveringNew")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_MIN_DIST_UNCVR, get_random_searcher_weights("MinDistToUncovered")));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_CVR_NEW, (double) this->coveredNew));
    dbl_map.insert(std::pair<std::string, double>(WEIGHT_FORK , (double) this->forkDisabled));

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

/**	Extracts a string representation in KQuery format of the current constraints, from the Execution state*/
const char * ExecutionState::get_constraints() {
    std::string Str;
    llvm::raw_string_ostream info(Str);
    ExprPPrinter::printConstraints(info, this->constraints);
    std::string res = info.str();

    return res.c_str();
}

/**
		Given an ExecutionState and a choosen WeightType in string format,
		Extract the weight given by WeightedRandomSearcher to this state with this weight type.
		Possible types are: Depth, InstCount, CPInstCount, QueryCost, CoveringNew, MinDistToUncovered
		*/
double ExecutionState::get_random_searcher_weights(const std::string type) {
    if (type == "Depth") {
        return this->weight;
    } else if (type == "InstCount") {
        uint64_t count = theStatisticManager->getIndexedValue(stats::instructions, this->pc->info->id);
        double inv = 1. / std::max((uint64_t) 1, count);
        return inv * inv;
    } else if (type == "CPInstCount") {
        StackFrame &sf = this->stack.back();
        uint64_t count = sf.callPathNode->statistics.getValue(stats::instructions);
        double inv = 1. / std::max((uint64_t) 1, count);
        return inv;
    } else if (type == "QueryCost") {
        return (this->queryCost < .1) ? 1. : 1./this->queryCost;
    } else if (type == "CoveringNew") {
        uint64_t md2u = computeMinDistToUncovered(this->pc, this->stack.back().minDistToUncoveredOnReturn);
        double invMD2U = 1. / (md2u ? md2u : 10000);
        double invCovNew = 0.;
        if (this->instsSinceCovNew)
            invCovNew = 1. / std::max(1, (int) this->instsSinceCovNew - 1000);
        return (invCovNew * invCovNew + invMD2U * invMD2U);
    } else if (type == "MinDistToUncovered") {
        uint64_t md2u = computeMinDistToUncovered(this->pc, this->stack.back().minDistToUncoveredOnReturn);
        double invMD2U = 1. / (md2u ? md2u : 10000);
        return invMD2U * invMD2U;
    }
    return -1.0;
}
