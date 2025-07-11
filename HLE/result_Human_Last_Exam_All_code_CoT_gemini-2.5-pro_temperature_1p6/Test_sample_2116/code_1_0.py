import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The derivation combines the properties of Pareto and LogSeries distributions,
    ultimately leading to the expression: E = pi / (2 * log(2)).
    """
    
    # Define the constants from the derived formula
    pi_val = math.pi
    # The natural logarithm of 2 is a key component
    log2_val = math.log(2)
    
    # The derived formula for the expected maximum magnitude
    expected_max_magnitude = pi_val / (2 * log2_val)

    # Output the explanation and the components of the final equation
    print("The final derived equation for the expected maximum earthquake magnitude E is:")
    print("E = pi / (2 * log(2))")
    print("\nHere are the values of each number in the equation:")
    print(f"pi = {pi_val}")
    print(f"2 = 2")
    print(f"log(2) = {log2_val}")

    print("\nPlugging these values into the equation:")
    print(f"E = {pi_val} / (2 * {log2_val})")
    
    print("\nThe calculated expected maximum earthquake magnitude is:")
    print(expected_max_magnitude)

solve_earthquake_magnitude()