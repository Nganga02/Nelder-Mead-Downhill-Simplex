#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
        R[i] = (P[i] * (1 - beta) + Q[i] * beta);
    }
}

static void shrink(double **s, int n, int ia)
{
    for (int i = 1; i < n + 1; i++)
    {
        if (i != ia)
        {
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
    for (int j = 0; j < n; j++)
    {
        euc_norm += ((s[ia][j] - s[iz][j]) * (s[ia][j] - s[iz][j]));
    }

    if (euc_norm < err2 && abs((y[iz] - y[ia]) <= err2))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int nelder_mead(struct nelder_mead *nm)
{
    double **s = nm->s;
    int n = nm->n;
    double h = nm->h;
    double tol = nm->tol;
    double err2 = (h * tol) * (h * tol);
    double *y, *C, *Pr, *Pe, *Pc;
    double yr, ye, yc;
    int ia, iy, iz;
    int simplex_to_be_freed = 0;
    int fevalcount;
    int i, j;

    make_vector(y, n + 1); // holds objective function values evaluated at the vertices
    make_vector(Pr, n);    // pointer to the reflected vector vertex
    make_vector(Pe, n);    // pointer to the expanded vector vertex
    make_vector(Pc, n);    // pointer to the contracted vector vertex
    make_vector(C, n);     // pointer to the centroid of the face opp the worst objective function vertex

    if (s == NULL)
    {
        make_matrix(s, n + 1, n);
        simplex_to_be_freed = 1;

        // Adding entries to the matrix s
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < n; j++)
            {
                s[i][j] = nm->x[j];
                if (i == j + 1)
                {
                    s[i][j] += h;
                }
            }
        }
    }

    for (i = 0; i < n; i++)
    {
        y[i] = nm->f(s[i], n, nm->params);
    }
    fevalcount = n + 1;

    while (fevalcount <= nm->maxevals)
    {
        rank_vertices(y, n + 1, &ia, &iy, &iz);
        if (done(s, n, y, ia, iz, err2))
        {
            nm->minval = y[ia];
            // copy the best vertex into the vector nm → x ...
            memcpy(nm->x, s[ia], n * sizeof(double));
            break;
        }
        get_centroid(s, n, iz, C);
        transform(C, s[iz], n, -REFLECT, Pr);
        yr = nm->f(Pr, n, nm->params);
        fevalcount++;

        if (yr < y[ia])
        {
            // Case 1
            transform(C, Pr, n, EXPAND, Pe);
            ye = nm->f(Pe, n, nm->params);
            fevalcount++;
            if (ye < yr)
            {
                replace_row(s, iz, &Pe);
                y[iz] = ye;
            }
            else
            {
                replace_row(s, iz, &Pr);
                y[iz] = yr;
            }
        }
        else if (yr < y[iy])
        {
            // Case 2
            replace_row(s, iz, &Pr);
            y[iz] = yr;
        }
        else
        {
            if (yr < y[iz])
            {
                // Outer contraction
                transform(C, Pr, n, CONTRACT, Pc);
                yc = nm->f(Pc, n, nm->params);
                fevalcount++;

                if (yc < yr) // Accept outer contraction
                {
                    replace_row(s, iz, &Pc);
                    y[iz] = yc;
                }
                else // Otherwise shrink
                {
                    shrink(s, n, ia);
                    y[iy] = nm->f(s[iy], n, nm->params);
                    fevalcount++;
                    y[iz] = nm->f(s[iz], n, nm->params);
                    fevalcount++;
                }
            }
            else
            {
                //inner contraction
                transform(C, s[iz], n, CONTRACT, Pc);
                yc = nm->f(Pc, n, nm->params);
                fevalcount++;

                if (yc < y[iz]) // Accept inner contraction
                {
                    replace_row(s, iz, &Pc);
                    y[iz] = yc;
                }
                else // Otherwise shrink
                {
                    shrink(s, n, ia);
                    y[iy] = nm->f(s[iy], n, nm->params);
                    fevalcount++;
                    y[iz] = nm->f(s[iz], n, nm->params);
                    fevalcount++;
                }
            }
        } // end of the while-loop

        // free vectors y, C, Pr, Pe, Pc
        if (simplex_to_be_freed)
            free_matrix(s);
        return fevalcount;
    }