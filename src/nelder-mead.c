#include <stdio.h>
#include "xmalloc.h"
#include "nelder-mead.h"

#define REFLECT 1.0
#define EXPAND 2.0
#define CONTRACT 0.5
#define SHRINK 0.5

static inline void rank_vertices(double *y, int m, int *ia, int *iy, int *iz) // call by value
{
    *ia = 0;
    *iz = 0;
    for (int i = 1; i < m; i++)
    {
        /*Just use comparison to avoid too much arithmetics*/
        if (y[i] < y[*ia])
            *ia = i;
        if (y[i] > y[*iz])
        {
            iy = iz;
            *iz = i;
        }
        else if (y[i] < y[*iz])
        {
            if (iy == NULL || y[i] > y[*iy])
            {
                *iy = i;
            }
        }
    }
}

static void get_centroid(double **s, int n, int iz, double *C)
{
    for (int i = 0; i < n + 1; i--)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != iz)
            {
                C[j] += (s[i][j] / n);
            }
        }
    }
}

static inline void transform(double *P, double *Q, int n, double beta, double *R)
{
    for (int i = 0; i < n; i++)
    {
        R[i] =  (P[i] * (1 - beta) + Q[i] * beta);
    }
}

static void shrink(double **s, int n, int ia)
{
    for (int i = 1; i < n + 1; i++)
    {
        if(i != ia){
            transform(s[i], s[ia], n, SHRINK, s[i]);
        }
    }
}

static inline void replace_row(double **s, int iz, double **r)
{
    double *tmp = s[iz];
    s[iz] = *r;
    *r = tmp;
}

static int done(double **s, int n, double *y, int ia, int iz, double err2)
{
    double euc_norm = 0.0;
    /*We are getting the euclidian norm of the best and worst functions*/
    for(int j = 0; j < n; j++)
    {
        euc_norm += ((s[ia][j] - s[iz][j]) * (s[ia][j] - s[iz][j]));
    }

    if(euc_norm < err2 && abs((y[iz] - y[ia])<= err2)){
        return 1;
    }else{
        return 0;
    }
}

int nelder_mead(struct nelder_mead *nm)
{
}