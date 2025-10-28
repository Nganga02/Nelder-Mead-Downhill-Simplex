#include <stdio.h>
#include <math.h>

#include "nelder-mead.h"

#define ROOT_TWO sqrt(2)
#define SQR(a) ((a) * (a))

// The main objective os to minimize the langragian #define ROOT_TWO sqrt(2)


struct NM_params
{
    void (*work_func)(double* x, int n, int *a, double c);
};

static void projection_func(double *x, int n, int *a, double c)
{
    double a_sqrd = 0;
    for (int i = 0; i < n; i++)
    {
        c -= a[i] * x[i];
        a_sqrd += (a[i] * a[i]);
    }

    // Avoiding division by zero
    if (a_sqrd == 0.0)
        return;
    for (int i = 0; i < n; i++)
    {
        x[i] += (c / a_sqrd) * a[i];
    }
}

static double constrain_func(double *x, int n, void *params)
{
    int a[] = {1, 1, 1};
    struct NM_params *proj_func = (struct NM_params *)params;
    void (*project_func)(double *, int, int *, double) = proj_func->work_func;
    project_func(x, n, a, 3.0);

    return SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
}

int main(void)
{
    double x[] = {1.0, 2.0, 3.0}; // the initial point
    int evalcount;

    struct NM_params ep = {
        .work_func = projection_func,
    };
    // Alert! C99-style initialization!
    struct nelder_mead NM = {
        .f = constrain_func, // the objective function
        .n = 3,              // the dimension of the space
        .s = NULL,           // delegate the construction of s
        .x = x,              // initial point / final point
        .h = 0.1,            // problem’s scale
        .tol = 1.0e-5,       // tolerance
        .maxevals = 1000,    // cap on function evaluations
        .params = &ep,       // no parameters
    };

    printf("Starting the evaluation");
    evalcount = nelder_mead(&NM);

    printf("\n Expected solution: min = 3 at (1, 1, 1)");
    if (evalcount > NM.maxevals)
    {
        printf("\nNo convergence after %d function evaluation",
               evalcount);
    }
    else
    {

       if (NM.n == 1)
    {
        printf("\nT\nComputed solution: min = %g at (%g)\n", NM.minval, x[0]);
    }
    else if (NM.n == 2)
    {
        printf("\nComputed solution: min = %g at  (%g, %g)", NM.minval, x[0], x[1]);
    }
    else if (NM.n == 3)
    {
        printf("\nComputed solution: min = %g at  (%g, %g, %g)", NM.minval, x[0], x[1], x[2]);
    }
    else
    {
        printf("\nComputed solution: min = %g at  (", NM.minval);
        for (int i = 0; i < NM.n; i++)
        {
            printf("%g", x[i]);
            if (i < NM.n - 1)
                printf(", ");
        }
        printf(")");
    }
        printf("\nconverged after %d function evaluations\n",
               evalcount);
    }
    return 0;
}