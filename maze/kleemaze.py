import random
import os

C_CODE = '''#include <stdio.h>
#include <stdlib.h>
#include <klee/klee.h>
#define H {}
#define W {}
#define ITERS ((W-1)/2)*((H-1)/2)
char maze[H][W] =  {};

void draw () {{
	int i, j;
	for (i = 0; i < H; i++){{
		for (j = 0; j < W; j++)
				  printf ("%c", maze[i][j]);
		printf ("\\n");
	}}
	printf ("\\n");
}}

int main (int argc, char *argv[]) {{
    int x, y;     //Player position
    int ox, oy;   //Old player position
    int i = 0;    //Iteration number
    char program[ITERS];
    x = 1;
    y = 1;
    maze[y][x] = 'X';
    klee_make_symbolic(program,ITERS,"program");
    while (i < ITERS) {{
        ox = x;    //Save old player position
        oy = y;
        //first advancment is through the door\wall:
        switch (program[i]) {{
            case 'w':
                y--;
                break;
            case 's':
                y++;
                break;
            case 'a':
                x--;
                break;
            case 'd':
                x++;
                break;
            default:
                printf("Wrong command!(only w,s,a,d accepted!)\\n");
                printf("You lose!\\n");
                exit (-1); //TODO: in linux should be exit
        }}
        //check if we hit a wall:
        if (maze[y][x] != ' ' ) {{
            printf("You lose\\n");
            exit (-2); 
        }} else {{ // if not we need to pass the wall, i.e. advance again:
            maze[y][x] = 'X';
            switch (program[i]) {{
                case 'w':
                    y--;
                    break;
                case 's':
                    y++;
                    break;
                case 'a':
                    x--;
                    break;
                case 'd':
                    x++;
                    break;
                default:
                    printf("Wrong command!(only w,s,a,d accepted!)\\n");
                    printf("You lose!\\n");
                    exit(-1); //TODO: in linux should be exit
            }}
        }}
        //finished advanding, check if we are there:
        if (maze[y][x] == '#') {{
            printf("You win!\\n");
            printf("Your solution \\n", program);
            klee_assert(0);  //Signal The solution!!
            exit(1);//TODO: in linux should be exit
        }}
        maze[y][x] = 'X';
        draw();          //draw it
        i++;
        sleep(1); //me wait to human
    }}
    printf("You lose!");
    return -1;
}}'''
ROW = 0
COL = 1

UNSET = -1
FREE = 0
BLOCK = 1

UP = 0
RIGHT = 1
DOWN = 2
LEFT = 3

VISITED = 1
TREASURE = 2

OUTPUT_FILE_NAME = "KleeMaze{}.c"

class Maze(object):

    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self.mat = [] #matrix of occupied / not occupied
        for i in range(rows):
            self.mat.append([0 for k in range(cols)])
        self.walls = []
        for i in range(rows):
            self.walls.append([[UNSET, UNSET, UNSET, UNSET] for k in range(cols)])
        self.make_wall_frame()
        self.treasure_marked = False

    def build(self, pos, new_pos):
        '''
        Build walls advancing from pos to new_pos
        barrier between rooms is freed,
        If new_pos is surrounded by visited cells, it is barricaded as well
        '''
        newrow, newcol = new_pos
        oldrow, oldcol = pos
        is_last = (len(self.get_free_neighbors(new_pos)) == 0)
        # print("bulding from ({x},{y}) to ({x2},{y2})".format(x=oldrow, y=oldcol, x2=newrow, y2=newcol))
        if newrow < oldrow: #up
            self.set_wall(pos, UP, FREE)
            self.set_wall(new_pos, DOWN, FREE)
            self.block_unset_walls(pos, (LEFT, DOWN, RIGHT))
            if is_last:
                self.block_unset_walls(new_pos, (UP, LEFT, RIGHT))
        elif newrow > oldrow: #down
            self.set_wall(pos, DOWN, FREE)
            self.set_wall(new_pos, UP, FREE)
            self.block_unset_walls(pos, (LEFT, UP, RIGHT))
            if is_last:
                self.block_unset_walls(new_pos, (DOWN, LEFT, RIGHT))
        elif newcol > oldcol: #right
            self.set_wall(pos, RIGHT, FREE)
            self.set_wall(new_pos, LEFT, FREE)
            self.block_unset_walls(pos, (UP, LEFT, DOWN))
            if is_last:
                self.block_unset_walls(new_pos, (RIGHT, UP, DOWN))
        elif newcol < oldcol: #left
            self.set_wall(pos, LEFT, FREE)
            self.set_wall(new_pos, RIGHT, FREE)
            self.block_unset_walls(pos, (UP, RIGHT, DOWN))
            if is_last:
                self.block_unset_walls(new_pos, (LEFT, UP, DOWN))
        else:
            raise IndexError("Bad building args")

    def make_wall_frame(self):
        for i in range(self.rows):
            self.set_wall((i, 0), LEFT, BLOCK)
            self.set_wall((i, self.cols-1), RIGHT, BLOCK)
        for j in range(self.cols):
            self.set_wall((0, j), UP, BLOCK)
            self.set_wall((self.rows-1, j), DOWN, BLOCK)

    def get_wall(self, pos, loc):
        '''
        Return the status of the wall at position and location
        '''
        return self.walls[pos[ROW]][pos[COL]][loc]

    def set_wall(self, pos, loc, val):
        self.walls[pos[ROW]][pos[COL]][loc] = val

    def block_unset_walls(self, pos, locs):
        '''
        put blocks on all walls on locs only if they are unset
        '''
        for loc in locs:
            if UNSET == self.get_wall(pos, loc):
                self.set_wall(pos, loc, BLOCK)

    def set(self, pos, val=VISITED):
        self.mat[pos[ROW]][pos[COL]] = val

    def get(self, pos):
        '''
        Safely get location. If out of bounds, returns None
        '''
        if 0 <= pos[ROW] < self.rows and 0 <= pos[COL]< self.cols:
            return self.mat[pos[ROW]][pos[COL]]
        else:
            return None

    def get_free_neighbors(self, pos):
        candidates = [(pos[ROW]+1, pos[COL]),
                      (pos[ROW]-1, pos[COL]),
                      (pos[ROW], pos[COL]+1),
                      (pos[ROW], pos[COL]-1)]
        return list(filter(lambda x: self.get(x) == FREE, candidates))

    def random_build(self, pos):
        '''
        Begins building randomly from pos until stuck (assumes starting pos is visited)
        '''
        opts = self.get_free_neighbors(pos)
        # print("Building from ({x},{y})".format(x=pos[ROW], y=pos[COL]))
        while len(opts) > 0:
            next_pos = random.choice(opts)
            self.set(next_pos) #set as visited
            self.build(pos, next_pos)
            pos = next_pos
            opts = self.get_free_neighbors(pos)
        if not self.treasure_marked:
            # print("Treasuer position: ({x},{y})".format(x=pos[ROW], y=pos[COL]))
            self.set(pos, TREASURE)
            self.treasure_marked = True
        return pos

    def get_random_start_pos(self):
        '''
        Sekects a random VISITED cell with free neigbor
        '''
        rows = list(range(self.rows))
        cols = list(range(self.cols))
        random.shuffle(rows)
        random.shuffle(cols)
        for i in rows:
            for j in cols:
                if self.get((i,j)) == VISITED and len(self.get_free_neighbors((i,j))) > 0:
                    return ((i,j))
        return None  # all full


    def get_c_style(self):
        rows = []
        row = 0
        while row < self.rows:
            #print top
            cur_row = ""
            for i in range(self.cols):
                w = self.get_wall((row, i), UP)
                if i == 0:
                    if w == BLOCK:
                        cur_row += '---'
                    else:
                        cur_row += '- -'
                else:
                    if w == BLOCK:
                        cur_row += '--'
                    else:
                        cur_row += ' -'
            rows.append(cur_row)
            #print middle
            cur_row = ""
            for j in range(self.cols):
                prev = None
                the_wall = self.walls[row][j]
                if (j ==0):
                    if (row == 0):
                        cur_row += "< "
                        if the_wall[RIGHT] == FREE:
                            cur_row += " "
                        else:
                            cur_row += "|"
                        continue
                    else:
                        cur_row+= "|"
                        if self.get((row,j)) == TREASURE:
                            cur_row += "#"
                        else:
                            cur_row += " "
                        if the_wall[RIGHT] == BLOCK:
                            cur_row += "|"
                        else:
                            cur_row += " "
                        continue
                if self.get((row,j)) == TREASURE:
                    cur_row += "#" #middle space
                else:
                    cur_row += " "
                if the_wall[RIGHT] == BLOCK:
                    cur_row += "|"
                else: cur_row += " "
            rows.append(cur_row)
            if row == self.rows-1:
                #print last bottom
                cur_row = "-"
                for i in range(self.cols):
                    cur_row += "--"
                rows.append(cur_row)
            row+=1
        return "{" + ",\n".join(map(lambda x : '\"' + x + '\"', rows)) + "}"


def get_random_maze(rows, cols):
    m = Maze(rows, cols)
    pos = (0, 0) #TODO: free the entrance
    m.set(pos)
    treasure = m.random_build(pos)
    pos = m.get_random_start_pos()
    while pos != None:
        m.random_build(pos)
        pos = m.get_random_start_pos()
    return m

for i in range(0, 100):
    m = get_random_maze(15,15)
    code = C_CODE.format(15*2+1, 15*2+1, m.get_c_style())
    with open(OUTPUT_FILE_NAME.format(i), 'w') as f:
        f.write(code)
    if (i ==0):
        print(os.path.abspath(OUTPUT_FILE_NAME.format(i)))
