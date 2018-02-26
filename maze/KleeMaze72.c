#include <stdio.h>
#include <stdlib.h>
#include <klee/klee.h>
#define H 31
#define W 31
#define ITERS ((W-1)/2)*((H-1)/2)
char maze[H][W] =  {"-------------------------------",
"<   |                   | |   |",
"--- --- ----------- --- - - - -",
"|#|   | |       |   |   |   | |",
"- --- --- --- - - --- ------- -",
"| |   |     | |   | |         |",
"- - --- ----- ----- -----------",
"| | |     |   |   |     |     |",
"- - ----- - --- - --- - - --- -",
"|       | |   | |   | |   |   |",
"- ------- --- - --- - ----- ---",
"| |       | |   |   |     |   |",
"- - - ----- ----- ----- - - - -",
"| | |     |     | |   | | | | |",
"- --- --- - ----- - - - - - ---",
"|   |   | |     |   |   | |   |",
"- - --- - ----- - --- ------- -",
"| | |   |   |   | |   |   |   |",
"- - - --- - - --- - ----- - ---",
"| |   | | | |     |       |   |",
"------- - - ----------------- -",
"|       | |   |       |     | |",
"- ------- --- - ----- - --- - -",
"|         |   | |   | |   | | |",
"--------------- - - - --- - - -",
"|   |           | | | |   |   |",
"- - - --- ------- - - - -------",
"| |   | |   |   | |   |   |   |",
"- ----- --- - - - ------- - - -",
"|         |   | |           | |",
"-------------------------------"};

void draw () {
	int i, j;
	for (i = 0; i < H; i++){
		for (j = 0; j < W; j++)
				  printf ("%c", maze[i][j]);
		printf ("\n");
	}
	printf ("\n");
}

int main (int argc, char *argv[]) {
    int x, y;     //Player position
    int ox, oy;   //Old player position
    int i = 0;    //Iteration number
    char program[ITERS];
    x = 1;
    y = 1;
    maze[y][x] = 'X';
    klee_make_symbolic(program,ITERS,"program");
    while (i < ITERS) {
        ox = x;    //Save old player position
        oy = y;
        //first advancment is through the door\wall:
        switch (program[i]) {
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
                printf("Wrong command!(only w,s,a,d accepted!)\n");
                printf("You lose!\n");
                exit (-1); //TODO: in linux should be exit
        }
        //check if we hit a wall:
        if (maze[y][x] != ' ' ) {
            printf("You lose\n");
            exit (-2); 
        } else { // if not we need to pass the wall, i.e. advance again:
            maze[y][x] = 'X';
            switch (program[i]) {
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
                    printf("Wrong command!(only w,s,a,d accepted!)\n");
                    printf("You lose!\n");
                    exit(-1); //TODO: in linux should be exit
            }
        }
        //finished advanding, check if we are there:
        if (maze[y][x] == '#') {
            printf("You win!\n");
            printf("Your solution \n", program);
            klee_assert(0);  //Signal The solution!!
            exit(1);//TODO: in linux should be exit
        }
        maze[y][x] = 'X';
        draw();          //draw it
        i++;
        sleep(1); //me wait to human
    }
    printf("You lose!");
    return -1;
}