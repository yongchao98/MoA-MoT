import math

def calculate_expected_max_magnitude():
    """
    This function calculates the expected maximum earthquake magnitude based on the
    derived formula: PI / (2 * log(2)).
    """
    # Define the mathematical constants from the final equation.
    # The expected maximum magnitude is given by pi / (2 * ln(2)).
    pi_val = math.pi
    denominator_constant = 2.0
    log_of_2 = math.log(2)

    # Calculate the result
    expected_magnitude = pi_val / (denominator_constant * log_of_2)

    # Output each number in the final equation as requested
    print("The final equation for the expected maximum magnitude is: PI / (2 * log(2))")
    print(f"Value of PI: {pi_val}")
    print(f"Value of the constant '2': {denominator_constant}")
    print(f"Value of natural logarithm of 2 (log(2)): {log_of_2}")
    
    # Print the equation with the numerical values
    final_equation_str = f"{pi_val} / ({denominator_constant} * {log_of_2})"
    print(f"The calculation is: {final_equation_str}")
    
    # Print the final answer
    print("\nResult:")
    print(f"The expected maximum earthquake magnitude is: {expected_magnitude}")

# Execute the function
calculate_expected_max_magnitude()