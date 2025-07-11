import numpy as np
from scipy.optimize import root_scalar

def solve_mle_cauchy():
    """
    Computes the MLE for the location parameter of a Cauchy distribution.
    """
    # The given sample S
    S = np.array([2.65, -6.91, 4.02, -5.97, 2.61])
    
    # Define the derivative of the log-likelihood w.r.t. theta.
    # We want to find the roots of this function.
    def log_likelihood_derivative(theta, data):
        return np.sum(2 * (data - theta) / (1 + (data - theta)**2))

    # The derivative may have multiple roots. We'll search for them in intervals
    # suggested by the data points themselves, as local extrema often occur near data.
    # Brackets are pairs [a,b] where the function has different signs at a and b.
    brackets = [
        (-7.0, -6.0),
        (2.6, 2.7),
        (4.0, 4.1)
    ]
    
    critical_points = []
    for bracket in brackets:
        try:
            sol = root_scalar(log_likelihood_derivative, args=(S,), bracket=bracket)
            critical_points.append(sol.root)
        except ValueError:
            # No root found in this bracket
            pass
            
    # Define the negative log-likelihood function (ignoring constants)
    # We want to find which critical point minimizes this function.
    def neg_log_likelihood(theta, data):
        return np.sum(np.log(1 + (data - theta)**2))

    # Evaluate the negative log-likelihood at each critical point
    nll_values = [neg_log_likelihood(p, S) for p in critical_points]
    
    # The MLE is the critical point that corresponds to the minimum NLL
    min_nll_value = min(nll_values)
    mle_theta = critical_points[nll_values.index(min_nll_value)]

    # --- Output Results ---
    
    print("The Maximum Likelihood Estimate (MLE) of theta is found by solving the equation where the derivative of the log-likelihood is zero:")
    print("d(log L)/d(theta) = sum( 2 * (xi - theta) / (1 + (xi - theta)^2) ) = 0")
    print("\nFor the sample S = {2.65, -6.91, 4.02, -5.97, 2.61}, we numerically find the MLE for theta.")
    print(f"\nThe MLE theta_hat is: {mle_theta:.7f}")
    
    print("\nSubstituting this value back into the equation term by term:")
    
    total_sum = 0
    equation_terms = []
    for x in S:
        term = 2 * (x - mle_theta) / (1 + (x - mle_theta)**2)
        equation_terms.append(f"({term:+.6f})")
        total_sum += term
        
    equation_str = " + ".join(equation_terms)
    print(f"  {equation_str} = {total_sum:.6e}")
    print("\nThe sum is close to zero, confirming this is a valid solution.")

    # Print the final answer rounded to one decimal place
    rounded_mle_theta = round(mle_theta, 1)
    print(f"\nThe Maximum Likelihood Estimate of theta is {mle_theta:.7f}")
    print(f"Rounded to one decimal place, the estimate is: {rounded_mle_theta}")
    
    return rounded_mle_theta

# Run the solver and get the final answer
final_answer = solve_mle_cauchy()

# The final answer in the required format
print(f"\n<<<{final_answer}>>>")
