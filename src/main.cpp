#include <iostream>
#include <vector>
#include <math.h>
#include "solver.hpp"
#include <GL/glut.h>
#include<chrono>

//using namespace std::chrono
//int N=100;
#define IX(i, j) ((i) + (N) * (j))

/*
  ----------------------------------------------------------------------
   MAIN
  ----------------------------------------------------------------------
*/

int main(int argc, char **argv)
{

  Domain domain;
  Case C;
  printf("\t ********************************************* \n");
  printf("\t *****-----Welcome to LBM Simulation-----***** \n");
  printf("\t ********************************************* \n");
  printf("\t Calculating Solution......... \n \n");
  printf("\t Iteration \t Error \n");
  auto start = std::chrono::steady_clock::now();
  C.Simulate(domain);
  auto end = std::chrono::steady_clock::now();
  
  auto t_elapsed=std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  printf("\t ******************************************* \n");
  printf("\t Simulation ends with time elapsed: %d ms. \n",t_elapsed);
  printf("\t ******************************************* \n");

  return 0;
}
