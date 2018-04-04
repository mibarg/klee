### for i in $(seq 0 $END); do /home/mibarg/source/klee/learn/mngr.sh KleeMaze$i.o 20 1000 100 ; done 
### command format: kleemaze.o max_time_per_run number_of_runs number_of_x0_per_run
# define maze to run
declare maze=~/source/klee/maze/$1
mazename="${1%.*}"
mkdir $mazename
cd $mazename
# run each opt separately
### Nelder-Mead"
declare opt="Nelder-Mead"
mkdir "WI-Nelder-Mead"
cp $maze "WI-Nelder-Mead"
cd "WI-Nelder-Mead"
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr $3 -ni 1 -wi -s -g "ICov(%)" -mm $opt
cd ..
### Powell
declare opt="Powell"
mkdir "WI-Powell"
cp $maze "WI-Powell"
cd "WI-Powell"
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr $3 -ni 1 -wi -s -g "ICov(%)" -mm $opt
cd ..
### Nelder-Mead"
declare opt="Nelder-Mead"
mkdir $opt
cp $maze $opt
cd $opt
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr $3 -ni $4 -s -g "ICov(%)" -mm $opt
cd ..
### Powell
declare opt="Powell"
mkdir $opt
cp $maze $opt
cd $opt
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr $3 -ni $4 -s -g "ICov(%)" -mm $opt
cd ..
### random
declare opt="random"
mkdir $opt
cp $maze $opt
cd $opt
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr $3 -ni $4 -s -g "ICov(%)" -mm $opt
cd ..
### BFS
mkdir bfs
cp $maze bfs
cd bfs
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr 1 -ni 1 -s -g "ICov(%)" -mm $opt -ks bfs
cd ..
### DFS
mkdir dfs
cp $maze dfs
cd dfs
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr 1 -ni 1 -s -g "ICov(%)" -mm $opt -ks dfs
cd ..
### random-path
mkdir random-path
cp $maze random-path
cd random-path
python3 ~/source/klee/learn/learn_weights.py -kd ~/klee_build_dir/bin -ki $1 -ka "--warnings-only-to-file -max-time=$2" -nr 1 -ni 1 -s -g "ICov(%)" -mm $opt -ks "random-path"
cd ..
###
cd ..
echo "done"
