//Kyle Pullen
//Newton and Secant Approximations

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
    //declaring all of my vars
    //first 5 are taking cmd line arguments
    double p0 = atof(argv[1]);
    double p1 = atof(argv[2]);
    double tol = atof(argv[3]);
    int max = atoi(argv[4]);
    int lim = atoi(argv[5]);
    //i'll use this stuff to build my polynomial
    int cubic_coefficient = atoi(argv[6]);
    int squared_coefficient = atoi(argv[7]);
    int first_coefficient = atoi(argv[8]);
    int constant = atoi(argv[9]);
    //incrementer and place holder
    int i = 1;
    double pnext = p0;
    double pprev = p0;
    double storage;

    //necessary print statements, just displaying what is passed as arguments in cmd line
    printf("Polynomial Root Finder by Kyle Pullen\n");
    printf("\nIpnut Parameters:\n");
    printf("\n\tp0 = %f\n", p0);
    printf("\tp1 = %f\n", p1);
    printf("\ttol = %f\n", tol);
    printf("\tmax = %i\n", max);
    printf("\nPolynomial is of order: %i\n\n", lim);
    printf("Terms of polynomial: %.1fx^%i + %.1fx^%i + %.1fx + %i\n\n",
            (float)cubic_coefficient, lim, (float)squared_coefficient, lim - 1, (float)first_coefficient, constant);

    //I built this function to evaluate my point at the given function
    double myFunction(double x)
    {
        return (cubic_coefficient* ( pow(x, lim) ) + squared_coefficient * (pow(x, lim - 1))
                 + first_coefficient * ( x ) + constant );
    }
    //I build this function to evaluate my point at the derivative of the given function
    double atDerivative(double x)
    {
        return (cubic_coefficient * (lim) * ( pow(x, lim - 1) ) + squared_coefficient * (lim - 1)* (pow(x, lim - 2))
                 + first_coefficient);
    }

    while(i < max)
    {
        //each time make sure I'm not trying to divide by 0
        if(atDerivative(pprev) != 0)
        {
            //this equation was given
            pnext = pprev - (myFunction(pprev) / atDerivative(pprev));
            //comparing to tol each time
            if(fabs(pnext - pprev) < tol)
            {
                //doing the newton approximation
                printf("\tp%i = %.15f\n", i, pnext);
                printf("\n\tSolution found after %i iterations: %.15f\n", i, pnext);
                break;

            }
            printf("\tp%i = %.15f\n", i, pnext);
            i++;
            pprev = pnext;
        }
        //if I hit this, derivative was zero
        else
        {
            printf("Unsuccessful, derivative is zero\n");
        }
        //if I reach this point, I've run out of allowed iterations
        if(i == max && abs(pnext- pprev) > tol)
        {
            printf("\nUnsuccessful, max number of iterations performed\n");
        }
    }


    printf("\nSecant Method:\n\n");
    //resetting my counter
    i = 2;
    pnext = p1;
    pprev = p0;

    while(i < max)
    {
        //preventing dividing by zero
        if((myFunction(pnext) - myFunction(pprev) == 0))
           {
               printf("Divide by zero error, exiting loop");
               break;
           }
        //making sure we're above tol
        if(fabs(pnext - pprev) >= tol)
        {
            //completing all of the computation and preparing values for the next iteration
            storage = (pnext - (myFunction(pnext) * (pnext - pprev)) / (myFunction(pnext) - myFunction(pprev)));
            pprev = pnext;
            pnext = storage;
            //printing what I found
            printf("\tp%i: %.15f\n", i, pnext);
            i++;
        }
        //When I reach this point, I've found my solution
        else if(fabs(pnext - pprev) < tol)
        {
            printf("\n\tSolution found after %i iterations: %.15f", i - 2, pnext);
            break;
        }

        if(i == max && abs(pnext- pprev) > tol)
        {
            printf("\nUnsuccessful, max number of iterations performed\n");
        }
    }

    return 0;
}

