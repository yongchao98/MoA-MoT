import math

def get_closed_formula(n):
    """
    This function provides the closed-form formula for P(n) for a given positive integer n.
    The derivation shows that P(n) = (2 * pi)^(n*(n-1)/4) / sqrt(n!).
    This script will print the formula with the specific numbers for the given n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the components of the formula
    exponent_numerator = n * (n - 1)
    exponent_denominator = 4
    
    # Simplify the exponent fraction for cleaner output
    common_divisor = math.gcd(exponent_numerator, exponent_denominator)
    exp_num_simple = exponent_numerator // common_divisor
    exp_den_simple = exponent_denominator // common_divisor

    if exp_den_simple == 1:
        exponent_str = f"{exp_num_simple}"
    else:
        exponent_str = f"{exp_num_simple}/{exp_den_simple}"

    n_factorial_val = math.factorial(n)

    # Print the final equation with the computed numbers
    print(f"For n = {n}, the closed-form formula for P(n) is:")
    print(f"P({n}) = (2 * pi)^({exponent_str}) / sqrt({n}!)")
    
    print("\nThe specific numbers in this equation are:")
    print(f"- The base of the power in the numerator is 2*pi.")
    print(f"- The exponent is {exponent_str}.")
    print(f"- The term inside the square root in the denominator, n!, is {n_factorial_val}.")

# Set the value of n for which to find the formula.
# You can change this value to any positive integer.
n = 5
get_closed_formula(n)