#include<iostream>
#include<vector>
#include <math.h>
//#include

struct Domain
{

    int imax;
    int jmax;
    int dim;
    int lattice;
    double dx;
    double dy;
    double nu;
    double omega;
    double Tlid;
    double F;
    double ulid;
    double Re;
    double err;

    // //int N;
    // float dt, diff, visc;
    // float force, source;
    // int dvel;

    // float *u, *v, *u_prev, *v_prev;
    // float *dens, *dens_prev;

    // int win_id;
    // int win_x, win_y;
    // int mouse_down[3];
    // int omx, omy, mx, my;

    std::vector<double> w;
    std::vector<double> cx;
    std::vector<double> cy;

    //std::vector<double> F;
    std::vector<std::vector<double>> RHO;
    std::vector<std::vector<double>> U;
    std::vector<std::vector<double>> VEL_OLD;
    std::vector<std::vector<double>> V;
    std::vector<std::vector<double>> VEL;
    std::vector<std::vector<double>> T;
    std::vector<std::vector<std::vector<double>>> f;
    std::vector<std::vector<std::vector<double>>> feq;

    Domain()
    {

        // win_x=512;
        // win_y=512;
        //N=100;
        imax = 51;
        jmax = 51;
        dim = 9;
        lattice;
        dx = 1;
        dy = 1;
        ulid=0.1;
        nu = 0.04;
        Tlid = 1.0;
        omega = 2 / (6 * nu + 1);
        F=0.01;
        
        // tau = 1 / omega;
        w = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
        cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
        cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};

        RHO.resize(imax, std::vector<double>(jmax, 5));
        U.resize(imax, std::vector<double>(jmax, 0));
        V.resize(imax, std::vector<double>(jmax, 0));
        VEL.resize(imax, std::vector<double>(jmax, 0));
        VEL_OLD.resize(imax, std::vector<double>(jmax, 0));
        T.resize(imax, std::vector<double>(jmax, 0));
        //F.resize(dim,0);
        f.resize(dim, std::vector<std::vector<double>>(imax, std::vector<double>(jmax, 0)));
        feq.resize(dim, std::vector<std::vector<double>>(imax, std::vector<double>(jmax, 0)));
        Re=1;
        //printf("Reynolds no: %f \n",ulid*imax/nu);
    }
};

/*
  ----------------------------------------------------------------------
    CLASS CASE
  ----------------------------------------------------------------------
*/

class Case
{

public:
  static  void Collision(Domain &domain);

  static  void Streaming(Domain &domain);

  static  void BoundaryConditions(Domain &domain);

  static  void Result(Domain &domain);

  static  void Simulate(Domain &domain);

  static double ErrCalc(Domain &domain);

    //void printmatrix(std::vector<std::vector<double>> &vec, int imax, int jmax);

private:
    int _q;
    int _nx;
    int _ny;
    Domain _domain;
};

