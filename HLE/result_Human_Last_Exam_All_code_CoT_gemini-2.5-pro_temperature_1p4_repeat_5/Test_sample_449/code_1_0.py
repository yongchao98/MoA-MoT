import math

def solve_probability():
    """
    Calculates the probability that the conditioned random walk never enters
    the set of the four neighbors of the origin.
    """
    # The starting point is (3000, 4000)
    norm_x0 = 5000.0

    # Euler-Mascheroni constant
    gamma_E = 0.5772156649015328

    # The formula for the probability is P = 1 - 2 / (ln(||x_0||) + gamma_E + ln(8))
    log_norm_x0 = math.log(norm_x0)
    log_8 = math.log(8)

    # Calculate the denominator of the fraction
    denominator = log_norm_x0 + gamma_E + log_8

    # Calculate the final probability
    prob = 1 - (2 / denominator)

    print("The probability P is calculated using the formula:")
    print("P = 1 - (2 / (ln(||x_0||) + gamma_E + ln(8)))")
    print("\nSubstituting the numerical values:")
    print(f"||x_0|| = {int(norm_x0)}")
    print(f"ln({int(norm_x0)}) = {log_norm_x0}")
    print(f"gamma_E = {gamma_E}")
    print(f"ln(8) = {log_8}")
    
    print("\nThe final equation with all numbers is:")
    print(f"P = 1 - (2 / ({log_norm_x0} + {gamma_E} + {log_8}))")
    print(f"P = 1 - (2 / {denominator})")
    print(f"P = 1 - {2 / denominator}")
    print(f"P = {prob}")

    # Format the result to two significant digits for the final answer
    prob_two_sig_digits = f"{prob:.2g}"
    print(f"\nThe approximate answer with two significant digits is {prob_two_sig_digits}.")

solve_probability()