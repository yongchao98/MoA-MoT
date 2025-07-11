import math

def solve_l(d):
    """
    Calculates the value of l(d).

    Args:
        d: An integer d >= 2.

    Returns:
        The exact value of l(d).
    """
    if d < 2:
        raise ValueError("d must be an integer greater than or equal to 2.")

    # The value of l(d) simplifies to (1/2) * ln((sqrt(d)-1)/(sqrt(d)+1))
    # which is equivalent to -arctanh(1/sqrt(d)).

    sqrt_d = math.sqrt(d)
    
    # We use the log form for calculation
    numerator = sqrt_d - 1
    denominator = sqrt_d + 1
    
    # Check for domain of log
    if numerator <= 0 or denominator <= 0:
        # This case happens for d=1, but d>=2 is given.
        # For d=2, numerator > 0.
        raise ValueError("Invalid value of d for log calculation.")

    val = 0.5 * math.log(numerator / denominator)
    
    return val

def print_final_equation(d, result):
    """
    Prints the final equation and its components as requested.
    """
    print(f"For d = {d}, the exact value is l({d}) = {result}")
    
    print("\nThe final equation is: l(d) = (1/2) * ln((sqrt(d) - 1) / (sqrt(d) + 1))")
    
    # As requested, printing each number in the final equation
    print("The numbers in the final equation are:")
    print("d =", d)
    print("Numerator constant: -1")
    print("Denominator constant: 1")
    print("Fraction multiplier: 1/2 (which consists of 1 and 2)")


# Example usage with d=4
d_value = 4
try:
    l_value = solve_l(d_value)
    print_final_equation(d_value, l_value)
except ValueError as e:
    print(e)

# Example usage with d=2
d_value = 2
try:
    l_value = solve_l(d_value)
    print_final_equation(d_value, l_value)
except ValueError as e:
    print(e)