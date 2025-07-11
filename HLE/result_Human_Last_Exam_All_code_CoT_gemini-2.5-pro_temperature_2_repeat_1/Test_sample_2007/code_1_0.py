import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate (MLE) for the parameter theta
    of a Cauchy distribution given a sample S.
    """
    # The given simple random sample
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # To find the MLE of theta, we need to maximize the log-likelihood function.
    # This is equivalent to minimizing the negative of the log-likelihood.
    # After removing constant terms, the function to minimize is:
    # g(theta) = sum(log(1 + (x_i - theta)^2)) for all x_i in S.
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # A good initial guess for the location parameter of a Cauchy distribution
    # is the sample median.
    initial_guess = np.median(S)

    # Use a numerical optimization method to find the value of theta
    # that minimizes the negative log-likelihood function.
    result = minimize(neg_log_likelihood, initial_guess, args=(S,))
    
    # The MLE is the value of theta found by the optimizer.
    mle_theta = result.x[0]

    # The problem asks for the estimate out to one decimal place.
    mle_theta_rounded = round(mle_theta, 1)

    print(f"The Maximum Likelihood Estimate (MLE) for theta is: {mle_theta_rounded}")
    print("-" * 50)
    
    # The MLE is the value of theta that solves the likelihood equation:
    # sum_{i=1 to 5} (x_i - theta) / (1 + (x_i - theta)^2) = 0
    # Here we demonstrate that our calculated MLE satisfies this equation.
    print(f"Verification using the unrounded MLE (theta = {mle_theta:.6f}):")
    print("Likelihood equation terms:")

    terms = (S - mle_theta) / (1 + (S - mle_theta)**2)
    total_sum = 0
    
    # Output each number/term in the final equation
    for x_i, term in zip(S, terms):
        print(f"  For x = {x_i:6.2f}, term = ({x_i:6.2f} - {mle_theta:.4f}) / (1 + ({x_i:6.2f} - {mle_theta:.4f})^2) = {term:10.7f}")
        total_sum += term

    print("-----------------------------------------------------------------------")
    print(f"Sum of terms = {total_sum:.10f}")
    print("The sum is approximately 0, which confirms our MLE is correct.")

if __name__ == '__main__':
    solve_mle()
    # The final answer required.
    print("\n<<<2.4>>>")