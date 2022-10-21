#include "solver.hpp"

/*
  ----------------------------------------------------------------------
   DOMAIN DECLARATION AND INITIALISATION
  ----------------------------------------------------------------------
*/

/*
  ----------------------------------------------------------------------
   PRINT MATRIX FUNCTION
  ----------------------------------------------------------------------
*/

void printmatrix(std::vector<std::vector<double>> &vec, int imax, int jmax)
{
    std::cout << " hello " << std::endl;
    for (int i = 0; i < vec.size(); i++)
    {

        for (int j = 0; j < vec[i].size(); j++)
        {
            std::cout << " " << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
  ----------------------------------------------------------------------
   INITIALIZE FROM CLASS CASE
  ----------------------------------------------------------------------
*/

/*
  ----------------------------------------------------------------------
   COLLISION FROM CLASS CASE
  ----------------------------------------------------------------------
*/

void Case::Collision(Domain &domain)
{
    //**************************COLLISION**************************//
    for (int i = 0; i < domain.imax; i++)
    {
        for (int j = 0; j < domain.jmax; j++)
        {
            for (int k = 0; k < domain.dim; k++)
            {
                double t1 = domain.cx[k] * domain.U[i][j] + domain.cy[k] * domain.V[i][j];
                double t2 = std::pow(domain.U[i][j], 2) + std::pow(domain.V[i][j], 2);

                domain.feq[k][i][j] = domain.w[k] * domain.RHO[i][j] * (1 + 3 * t1 + 4.5 * std::pow(t1, 2) - 1.5 * t2);
                domain.f[k][i][j] = domain.f[k][i][j] * (1 - domain.omega) + domain.omega * domain.feq[k][i][j]; // + domain.F;
            }
        }
    }
}

/*
  ----------------------------------------------------------------------
   STREAMING FROM CLASS CASE
  ----------------------------------------------------------------------
*/

void Case::Streaming(Domain &domain)
{

    //**************************STREAMING**************************//

    // for (int i = 0; i < domain.imax; i++) // j=1:ny
    // {
    //     for (int j = 0; j < domain.jmax; j++) // i=nx:-1:2
    //     {

    //         //        yn=mod(j,ny)+1;
    //         //        xe=mod(i,nx)+1;
    //         //        ys=ny-mod(ny+1-j,ny);
    //         //        xw=nx-mod(nx+1-i,ny);

    //         double yn = (j % domain.jmax-1) + 1;
    //         double xe = (i % domain.imax-1) + 1;

    //         double ys = domain.jmax-1 - (domain.jmax-1 + 1 - j) % domain.jmax-1;
    //         double xw = domain.imax-1 - (domain.imax-1 + 1 - i) % domain.imax-1;

    //         domain.feq[1][xe][j] = domain.f[1][i][j];
    //         domain.feq[2][i][yn] = domain.f[2][i][j];
    //         domain.feq[3][xw][j] = domain.f[3][i][j];
    //         domain.feq[4][i][ys] = domain.f[4][i][j];
    //         domain.feq[5][xe][yn] = domain.f[5][i][j];
    //         domain.feq[6][xw][yn] = domain.f[6][i][j];
    //         domain.feq[7][xw][ys] = domain.f[7][i][j];
    //         domain.feq[8][xe][ys] = domain.f[8][i][j];
    //     }
    // }

    for (int j = 0; j < domain.jmax; j++)         // j=1:ny
        for (int i = domain.imax - 1; i > 0; i--) // i=nx:-1:2
            domain.f[1][i][j] = domain.f[1][i - 1][j];

    for (int j = domain.jmax - 1; j > 0; j--) // ny:-1:2
        for (int i = 0; i < domain.imax; i++) // 1:nx
            domain.f[2][i][j] = domain.f[2][i][j - 1];

    for (int j = 0; j < domain.jmax; j++)         // j=1:ny
        for (int i = 0; i < domain.imax - 1; i++) // i=1:nx-1
            domain.f[3][i][j] = domain.f[3][i + 1][j];

    for (int j = 0; j < domain.jmax - 1; j++) // j=1:ny-1
        for (int i = 0; i < domain.imax; i++) // i=1:nx
            domain.f[4][i][j] = domain.f[4][i][j + 1];

    for (int j = domain.jmax - 1; j > 0; j--)     // ny:-1:2
        for (int i = domain.imax - 1; i > 0; i--) // nx:-1:2
            domain.f[5][i][j] = domain.f[5][i - 1][j - 1];

    for (int j = domain.jmax - 1; j > 0; j--)     // j=ny:-1:2
        for (int i = 0; i < domain.imax - 1; i++) // i=1:nx-1
            domain.f[6][i][j] = domain.f[6][i + 1][j - 1];

    for (int j = 0; j < domain.jmax - 1; j++)     // j=1:ny-1
        for (int i = 0; i < domain.imax - 1; i++) // i=1:nx-1
            domain.f[7][i][j] = domain.f[7][i + 1][j + 1];

    for (int j = 0; j < domain.jmax - 1; j++)     // j=1:ny-1
        for (int i = domain.imax - 1; i > 0; i--) // i=nx:-1:2
            domain.f[8][i][j] = domain.f[8][i - 1][j + 1];
}

/*
  ----------------------------------------------------------------------
   BOUNDARY CONDITIONS FROM CLASS CASE
  ----------------------------------------------------------------------
*/

void Case::BoundaryConditions(Domain &domain)
{

    //  Bottom BC
    for (int i = 0; i < domain.imax; i++)
    {
        domain.f[6][i][0] = domain.f[8][i][0];
        domain.f[2][i][0] = domain.f[4][i][0];
        domain.f[5][i][0] = domain.f[7][i][0];
    }

    //  Right BC
    for (int j = 0; j < domain.jmax; j++) // j=1:ny
    {
        domain.f[3][domain.imax - 1][j] = domain.f[1][domain.imax - 1][j];
        domain.f[7][domain.imax - 1][j] = domain.f[5][domain.imax - 1][j];
        domain.f[6][domain.imax - 1][j] = domain.f[8][domain.imax - 1][j];
    }

    //  Left BC
    for (int j = 0; j < domain.jmax; j++) // j=1:ny
    {
        domain.f[5][0][j] = domain.f[7][0][j];
        domain.f[8][0][j] = domain.f[6][0][j];
        domain.f[1][0][j] = domain.f[3][0][j];
    }

    // Top BC

    for (int i = 0; i < domain.imax; i++)
    {
        double rhon = domain.f[0][i][domain.jmax - 1] + domain.f[1][i][domain.jmax - 1] + domain.f[3][i][domain.jmax - 1] + 2 * (domain.f[2][i][domain.jmax - 1] + domain.f[5][i][domain.jmax - 1] + domain.f[6][i][domain.jmax - 1]);

        domain.f[7][i][domain.jmax - 1] = domain.f[5][i][domain.jmax - 1] + 0.5 * (domain.f[1][i][domain.jmax - 1] - domain.f[3][i][domain.jmax - 1]) - (0.5) * rhon * domain.ulid;
        domain.f[8][i][domain.jmax - 1] = domain.f[6][i][domain.jmax - 1] + 0.5 * (domain.f[3][i][domain.jmax - 1] - domain.f[1][i][domain.jmax - 1]) + (0.5) * rhon * domain.ulid;
        domain.f[4][i][domain.jmax - 1] = domain.f[2][i][domain.jmax - 1];

        // domain.f[6][i][domain.jmax - 1] =  -domain.f[8][i][domain.jmax - 1];
        // domain.f[4][i][domain.jmax - 1] =  -domain.f[2][i][domain.jmax - 1];
        // domain.f[7][i][domain.jmax - 1] =  -domain.f[5][i][domain.jmax - 1];
    }
}

/*
  ----------------------------------------------------------------------
   RESULT FROM CLASS CASE
  ----------------------------------------------------------------------
*/

void Case::Result(Domain &domain)
{
    for (int i = 0; i < domain.imax; i++)
    {
        for (int j = 0; j < domain.jmax; j++)
        {
            double sum = 0;
            for (int k = 0; k < domain.dim; k++)
            {
                sum += domain.f[k][i][j];
            }
            domain.RHO[i][j] = sum;
        }
    }

    for (int i = 0; i < domain.imax; i++)
    {
        for (int j = 0; j < domain.jmax; j++)
        {
            double sum2 = 0, sum3 = 0;
            for (int k = 0; k < domain.dim; k++)
            {
                sum2 += domain.f[k][i][j] * domain.cx[k];
                sum3 += domain.f[k][i][j] * domain.cy[k];
            }
            domain.U[i][j] = sum2 / domain.RHO[i][j];
            domain.V[i][j] = sum3 / domain.RHO[i][j];
        }
    }

    for (int i = 0; i < domain.imax; i++)
    {
        for (int j = 0; j < domain.jmax; j++)
        {

            domain.VEL[i][j] = std::pow(std::pow(domain.U[i][j], 2) + std::pow(domain.V[i][j], 2), 0.5);
        }
    }
    //error=ErrCalc(domain);
}

/*
  ----------------------------------------------------------------------
   SIMULATE FROM CLASS CASE
  ----------------------------------------------------------------------
*/

void Case::Simulate(Domain &domain)
{
    int iteration = 0;
    double error = 0;
    int maxiteration = 10000;

    // for (int i = 0; i < domain.imax; i++)
    //     domain.U[i][domain.jmax - 1] = domain.ulid;

    while (iteration < maxiteration)
    {
        domain.VEL_OLD = domain.VEL;
        Collision(domain);
        Streaming(domain);
        BoundaryConditions(domain);
        Result(domain);
        error = ErrCalc(domain);

        iteration++;

        if (iteration % 50 == 0)
            printf("\t   %d \t   %d \n", iteration, error);
    }
    // printmatrix(domain.VEL, domain.imax, domain.jmax);
}

double Case::ErrCalc(Domain &domain) // errcalculate(Domain &domain)
{
    int sum = 0;
    double error=0;
    for (int i = 1; i < domain.imax - 1; i++)
    {
        for (int j = 1; j < domain.jmax - 1; j++)
        {
            error=std::fabs(domain.VEL[i][j] - domain.VEL_OLD[i][j]);
            //sum = sum + std::pow((domain.VEL[i][j] - domain.VEL_OLD[i][j]), 2);
        }
    }
    
    //error = std::pow(sum, 0.5) / (domain.imax * domain.jmax);
    return error;
}