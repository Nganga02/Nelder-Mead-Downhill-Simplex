#include <stdio.h>
#include <math.h>
#include "nelder-mead.h"

#define ROOT_TWO sqrt(2)


static double w_func(double x);

struct en_params{
    double (*work_func)(double x);
    const double *load_vector;
};

static double w_func(double x)
{
    double x_sqrd = x*x;
    return ((x_sqrd*x_sqrd)/24) + (1/(12*x_sqrd)) - 0.125;
}

static double en_func(double *x, int n, void *params)
{
    struct en_params *p = (struct en_params *) params;
    double (*work_func)(double ) = p->work_func;
    const double *lv = p->load_vector;

    double term_1 = sqrt((x[0] * x[0]) + (x[1] * x[1])) / ROOT_TWO;
    double term_2 = sqrt((x[0] * x[0]) + (x[1] * x[1]));
    return (ROOT_TWO * work_func(term_1)) + work_func(term_2) - (lv[0] * (x[0] - 1)) - (lv[1] * (x[1] - 1));
}

int main(void)
{
    const double load_vector[] = {0, -0.3};
    double x[] = {-1.0, 1.0}; // the initial point
    int evalcount;

     struct en_params ep = {
        .work_func = w_func,
        .load_vector = load_vector
    };
    // Alert! C99-style initialization!
    struct nelder_mead NM = {
        .f = en_func,     // the objective function
        .n = 2,           // the dimension of the space
        .s = NULL,        // delegate the construction of s
        .x = x,           // initial point / final point
        .h = 0.1,         // problem’s scale
        .tol = 1.0e-4,    // tolerance
        .maxevals = 1000, // cap on function evaluations
        .params = &ep,   // no parameters
    };

    printf("Starting the evaluation");
    evalcount = nelder_mead(&NM);

    printf("\n Expected solution: min = -0.6395281605 at (0.8319642234, -1.2505278260)");
    if (evalcount > NM.maxevals)
    {
        printf("\nNo convergence after %d function evaluation",
               evalcount);
    }
    else
    {

        printf("\nComputed solution: min = %g at (%g, %g)\n",
               NM.minval, x[0], x[1]);
        printf("\nconverged after %d function evaluations",
               evalcount);
    }
    return 0;
}
