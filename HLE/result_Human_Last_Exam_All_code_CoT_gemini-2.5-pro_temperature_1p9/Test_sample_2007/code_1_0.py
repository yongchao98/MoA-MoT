import numpy as np
from scipy.optimize import minimize

def compute_mle():
    """
    This script computes the Maximum Likelihood Estimate (MLE) of the parameter theta
    for a Cauchy distribution given a sample S.
    """
    # The sample data provided
    S = [2.65, -6.91, 4.02, -5.97, 2.61]

    # The MLE for theta is found by maximizing the log-likelihood function.
    # This is equivalent to minimizing the negative log-likelihood.
    # For a Cauchy distribution, this simplifies to minimizing the sum of log(1 + (x_i - theta)^2).
    
    # Define the function to be minimized
    def neg_log_likelihood(theta, data):
        """
        Calculates the part of the negative log-likelihood to be minimized.
        """
        return np.sum(np.log(1 + (np.array(data) - theta)**2))

    # A robust initial guess for the location parameter of a Cauchy distribution
    # is the sample median.
    initial_guess = np.median(S)

    # Perform the minimization to find the MLE for theta
    result = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        args=(S,),
        method='Nelder-Mead'
    )

    # The MLE is the value of theta found by the optimizer
    mle_theta = result.x[0]

    # The MLE solves the equation where the derivative of the log-likelihood is zero:
    # Sum_{i} (x_i - theta) / (1 + (x_i - theta)^2) = 0
    # The following code prints this equation with the computed values.
    print("The maximum likelihood estimate of theta solves the equation where the derivative of the log-likelihood is zero.")
    print("The equation is: Sum[(x_i - theta) / (1 + (x_i - theta)^2)] = 0")
    print("Substituting the sample data and the estimated theta value:")
    
    equation_terms = []
    sum_of_terms = 0.0
    for x_i in S:
        term_numerator = x_i - mle_theta
        term_denominator = 1 + (x_i - mle_theta)**2
        term_value = term_numerator / term_denominator
        sum_of_terms += term_value
        equation_terms.append(f"({x_i:.2f} - {mle_theta:.4f}) / (1 + ({x_i:.2f} - {mle_theta:.4f})^2)")

    full_equation_str = " + ".join(equation_terms)
    
    print(f"\n{full_equation_str} = {sum_of_terms:.4e}\n")

    # The computed MLE value, rounded to one decimal place as requested.
    mle_theta_rounded = round(mle_theta, 1)

    print(f"The unrounded Maximum Likelihood Estimate for theta is: {mle_theta:.4f}")
    print(f"The Maximum Likelihood Estimate for theta rounded to one decimal place is: {mle_theta_rounded}")

compute_mle()
<<<2.6>>>