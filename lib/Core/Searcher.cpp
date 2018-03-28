//===-- Searcher.cpp ------------------------------------------------------===//
//
//                     The KLEE Symbolic Virtual Machine
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "Searcher.h"
#include <dirent.h>
#include "CoreStats.h"
#include "Executor.h"
#include "PTree.h"
#include <unistd.h>
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

#include "llvm/Support/Errno.h"
#include "klee/Constraints.h"
//#include "klee/util/ExprPPrinter.h" \\TODO
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>
#include <climits>
#include <sys/socket.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sstream>


using namespace klee;
using namespace llvm;
using namespace std;
using namespace std::chrono;

namespace {
  cl::opt<bool>
  DebugLogMerge("debug-log-merge");
}

namespace klee {
  extern RNG theRNG;
}

//-----------------STATIC MEMBERS definition:
std::string SMARTSearcher::outfilePath ;
std::string SMARTWeightedSearcher::weightsFile = "weights.json";
std::string SMARTWeightedSearcher::outfilePath;
unsigned int NeuralNetSearcher::destPort = 2805;
std::mutex NeuralNetSearcher::m;




Searcher::~Searcher() {
}

//-----------------SMARTSearcher:-----------------//
ExecutionState &SMARTSearcher::selectState() {
  return *states[theRNG.getInt32()%states.size()];
}
void SMARTSearcher::setOutputFileName(){
  char cwd[1024];
  std::string prefix = "klee-out-";
  int counter = 0;
  bool shouldUseTime = false;unsigned char isFolder = 0x4;
    //Failed getting pwd, write to a file with a name of time in miliseconds:
  if (!(getcwd(cwd, sizeof(cwd)) != NULL)) {     shouldUseTime = true;  }
  //write to states#.log with # current run according to the latest klee-out-#/
  else {
      std::string pwd(cwd);
      DIR *dir;
      struct dirent *ent;
      dir = opendir(pwd.c_str());
      if (dir != NULL) {
          while ((ent = readdir(dir)) != NULL) {
              auto res = std::mismatch(prefix.begin(), prefix.end(), std::string(ent->d_name).begin());
              if (ent->d_type == isFolder && res.first == prefix.end()) {
                  counter++;
              }
          }
          closedir(dir);
          SMARTSearcher::outfilePath = pwd + "/klee-last/" + "out.log";
      } else {    /* could not open directory */
          shouldUseTime = false;
      }
  }
  if (shouldUseTime){
    milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    SMARTSearcher::outfilePath = std::to_string(ms.count())+".log";
  }

}

void SMARTSearcher::update(ExecutionState *current,const std::vector<ExecutionState *> &addedStates,const std::vector<ExecutionState *> &removedStates) {

  states.insert(states.end(),addedStates.begin(), addedStates.end());
  //Write the new states to the file:
  std::ofstream file;

  for (auto state : addedStates){
      std::map<std::string, double> dbl_map;
      state->extractStateMaps(dbl_map, false);
      if (firstInsert == true){
        SMARTSearcher::setOutputFileName();
        std::cout<<"Writing states to file: "<< outfilePath<< std::endl;
        file.open(outfilePath, std::fstream::app);
          for(auto tup:dbl_map) {
              file<<tup.first << ",";

          }
          firstInsert  = false;
        file<<endl;
      } else{
        file.open(outfilePath, std::fstream::app);
      }
    for(auto tup:dbl_map) {
      file<<tup.second << ",";

    }
    file<<std::endl;
  }
  file.close();
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


void SMARTWeightedSearcher::setOutputFileName(){
    char cwd[1024];
    std::string prefix = "klee-out-";
    int counter = 0;
    bool shouldUseTime = false;unsigned char isFolder = 0x4;
    //Failed getting pwd, write to a file with a name of time in miliseconds:
    if (!(getcwd(cwd, sizeof(cwd)) != NULL)) {     shouldUseTime = true;  }
        //write to states#.log with # current run according to the latest klee-out-#/
    else {
        std::string pwd(cwd);
        DIR *dir;
        struct dirent *ent;
        dir = opendir(pwd.c_str());
        if (dir != NULL) {
            while ((ent = readdir(dir)) != NULL) {
                auto res = std::mismatch(prefix.begin(), prefix.end(), std::string(ent->d_name).begin());
                if (ent->d_type == isFolder && res.first == prefix.end()) {
                    counter++;
                }
            }
            closedir(dir);
            SMARTSearcher::outfilePath = pwd + "/" + std::to_string(counter - 1) + ".log";
        } else {    /* could not open directory */
            shouldUseTime = false;
        }
    }
    if (shouldUseTime){
        milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
        SMARTSearcher::outfilePath = std::to_string(ms.count())+".log";
    }

}

//-----------------SMARTWeightedSearcher:-----------------//
SMARTWeightedSearcher::SMARTWeightedSearcher(): states(new DiscretePDF<ExecutionState*>()){readWeights();}

SMARTWeightedSearcher::~SMARTWeightedSearcher(){ delete states;};

ExecutionState& SMARTWeightedSearcher::selectState(){return *states->choose(theRNG.getDoubleL());};

bool SMARTWeightedSearcher::empty() { return states->empty(); };

void SMARTWeightedSearcher::readWeights(){
  std::ifstream wFile (weightsFile);
  std::string line;
  bool finished = false;
  std::size_t lastChar = 0;
  if(wFile.is_open()){
    std::getline(wFile, line);
    wFile.close();
    //parse the json:
    while(!finished){
      //first find key:
      std::size_t firstApos = line.find("\"", lastChar);
      std::size_t secondApos = line.find("\"", firstApos+1);
      std::size_t startDbl = line.find(" ", secondApos + 1);
      std::size_t endDbl = line.find(",", startDbl + 1);
      std::string key = line.substr(firstApos+1, secondApos-firstApos -1 ); //TODO is this exactly?
      if (endDbl == std::string::npos){
        endDbl = line.find("}", startDbl);
        finished = true;
      }
      double value = std::stod(line.substr(startDbl + 1, endDbl - startDbl - 1));
      weights.emplace(key, value);

        std::cout << key << " : " << value << std::endl;
        lastChar = endDbl;
    }
      std::cout<< "added weight from file"<< std::endl;
  } else{
      weights.emplace(INSTR_PC_OPCODE,1);
      weights.emplace(INSTR_PC_OPERANDS,1);
    weights.emplace(INSTR_PC_USE_EMPTY,1);
    weights.emplace(INSTR_PPC_OPERANDS ,1);
    weights.emplace(INSTR_PPC_OPCODE ,1);
    weights.emplace(INSTR_PPC_USE_EMPTY,1);
    weights.emplace(SAT_NUM_ANDS,1);
    weights.emplace(SAT_NUM_ASSERTED,1);
    weights.emplace(SAT_NUM_EQUAL ,1);
    weights.emplace(SAT_NUM_XOR,1);
    weights.emplace(SAT_NUM_BOOL,1);
    weights.emplace(SAT_NUM_THEORY,1);
    weights.emplace(SAT_NUM_CLAUSES_VARS,1);
    weights.emplace(SAT_NUM_VARS_CLAUSES,1);
    weights.emplace(SAT_NUM_SUB,1);
    weights.emplace(SAT_NUM_MUL,1);
    weights.emplace(SAT_NUM_UDIV,1);
    weights.emplace(SAT_NUM_SDIV,1);
    weights.emplace(SAT_NUM_UREM,1);
    weights.emplace(SAT_NUM_OR,1);
    weights.emplace(SAT_NUM_SHL,1);
    weights.emplace(SAT_NUM_LSHR,1);
    weights.emplace(SAT_NUM_ASHR,1);
    weights.emplace(SAT_NUM_EXTRACT,1);
    weights.emplace(SAT_NUM_CONCAT,1);
    weights.emplace(SAT_NUM_ZEXT,1);
    weights.emplace(SAT_NUM_SEXT,1);
    weights.emplace(SAT_NUM_MSB,1);
    weights.emplace(SAT_NUM_MSB,1);
    weights.emplace(SAT_CLS_STD,1);
    weights.emplace(SAT_CLS_AVG,1);
    weights.emplace(SAT_CLS_MAX,1);
    weights.emplace(SAT_NUM_CAPS,1);
    weights.emplace(WEIGHT_DEPTH,1);
    weights.emplace(WEIGHT_INST_CNT,1);
    weights.emplace(WEIGHT_CPI_CNT,1);
    weights.emplace(WEIGHT_QUERY_CST,1);
    weights.emplace(WEIGHT_COVER_NEW,1);
    weights.emplace(WEIGHT_MIN_DIST_UNCVR,1);
    weights.emplace(WEIGHT_CVR_NEW,1);
    weights.emplace(WEIGHT_FORK,1);
    std::cout << "used default weights" << std::endl;
  }
}

double SMARTWeightedSearcher::getWeight(ExecutionState* state){
  std::map<std::string, double> dblmap;
  double score = 0;
  state->extractStateMaps(dblmap, 0);
  for (auto tup : dblmap){
    score += tup.second*weights[tup.first];
  }
  return score;
}

void SMARTWeightedSearcher::update(ExecutionState *current,const std::vector<ExecutionState *> &addedStates,const std::vector<ExecutionState *> &removedStates){
    if (current && std::find(removedStates.begin(), removedStates.end(), current) == removedStates.end()){
        states->update(current, getWeight(current));
    }
    //add states to data structure:
    std::ofstream file;
    for (std::vector<ExecutionState *>::const_iterator it = addedStates.begin(),ie = addedStates.end();it != ie; ++it) {
    ExecutionState *state = *it;
    states->insert(state, getWeight(state));
      //write to file:
      std::map<std::string, double> dbl_map;
      state->extractStateMaps(dbl_map, false);
      //Write headers:
      if (firstInsert == true){
          SMARTWeightedSearcher::setOutputFileName();
          std::cout<<"Writing states to file: "<< outfilePath<< std::endl;
          file.open(outfilePath, std::fstream::app);
          for(auto tup:dbl_map) {file<<tup.first << ",";}
          firstInsert  = false;
          file<<endl;
      } else{ file.open(outfilePath, std::fstream::app);}
      for(auto tup:dbl_map) {
          file<<tup.second << ",";

      }
      file<<std::endl;
  }
    file.close();
    //Remover states:
  for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
               ie = removedStates.end();
       it != ie; ++it) {
    states->remove(*it);
  }
}

//-----------------NeuralNetSearcher:-----------------//
NeuralNetSearcher::NeuralNetSearcher(): states(new DiscretePDF<ExecutionState*>()){openSocket();}

NeuralNetSearcher::~NeuralNetSearcher(){ delete states; close(socketFd);};

ExecutionState& NeuralNetSearcher::selectState(){return *states->choose(theRNG.getDoubleL());};

bool NeuralNetSearcher::empty() { return states->empty(); };

void NeuralNetSearcher::openSocket(){
    struct sockaddr_in serv_addr;
    struct timeval timeout;
    timeout.tv_sec = 1;
    timeout.tv_usec = 0;
    if ((socketFd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    {
        std::cout<<"Error connecting to neural network"<<std::endl;
        return;
    }
    memset(&serv_addr, '0', sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(NeuralNetSearcher::destPort);
    if(inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr)<=0)
    {
        std::cout<<"Error connecting to neural network"<<std::endl;
        return;
    }
    if (connect(socketFd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    {
        std::cout<<"Error connecting to neural network"<<std::endl;
        std::cout<<"Error connecting to neural network"<<std::endl;
        return;
    }
    if (setsockopt (socketFd, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout)) < 0)
        return;
    if (setsockopt (socketFd, SOL_SOCKET, SO_SNDTIMEO, (char *)&timeout, sizeof(timeout)) < 0)
        return;
    socketOpen = true;
    std::cout <<"opened socket"<< std::endl;
}

void NeuralNetSearcher::update(ExecutionState *current, const std::vector<ExecutionState *> &addedStates,const std::vector<ExecutionState *> &removedStates) {
    if (current && std::find(removedStates.begin(), removedStates.end(), current) == removedStates.end()){
        states->update(current, getWeight(current));
    }
    //add states to data structure:
    for (std::vector<ExecutionState *>::const_iterator it = addedStates.begin(),ie = addedStates.end();it != ie; ++it) {
        ExecutionState *state = *it;
        states->insert(state, getWeight(state));
    }
    //Remover states:
    for (std::vector<ExecutionState *>::const_iterator it = removedStates.begin(),
                 ie = removedStates.end();
         it != ie; ++it) {
        states->remove(*it);
    }
}

double NeuralNetSearcher::getWeight(ExecutionState * state) {
    if (!socketOpen){
        return 1.1;
    }
    std::stringstream json;
    std::map<std::string, double> dbl_map;
    state-> extractStateMaps(dbl_map, false);
    bool first = true;
    json << "{";
    for(auto tup:dbl_map) {
        if (!first){
            json <<", ";
        }
        json << "\"" << tup.first << "\"" << ": " << tup.second;
        first = false;
    }
    json << "}";
    //start mutex lock
//    std::unique_lock<std::mutex> lk(m);
    int res = write(socketFd, json.str().c_str(),json.str().size());
    int total = res;
    while (total< json.str().size() && res >0 ){
        res = write(socketFd, &(json.str().c_str()[total]),json.str().size() - total);
        total += res;
    }
    if (total<0){
  //      lk.unlock();
        socketOpen = false;
        return 1.1;
    }
    char buff[33];
    res = read(socketFd, buff, 32);
    total = res;
    while (res >0 && total < 32){
        res = read(socketFd, &(buff[total]), 32 - total);
        total += res;
    }
    if (total < 0){
        socketOpen = false;
    }
    if (total<=0){
    //    lk.unlock();
        return 1.1;
    }
    //lk.unlock();
    //end lock
    buff[total] = '\0';
    return std::stod(std::string(buff));
}

//-----------------End of compilation project-----------------//

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
             ie =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                states.end(); it != ie; ++it) {
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
	int is_debug = 1;
	states.front()->extractStateMaps(dbl_map,is_debug);
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
