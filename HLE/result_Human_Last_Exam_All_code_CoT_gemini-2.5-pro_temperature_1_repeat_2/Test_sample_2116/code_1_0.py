import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the problem's distributions.
    
    The final derived formula for the expected maximum magnitude E[M] is:
    E[M] = pi / (2 * log(2))
    
    This function calculates this value and prints the components of the formula
    and the final result.
    """
    
    # Define the constants from the formula
    pi_val = math.pi
    two_val = 2.0
    log2_val = math.log(2)
    
    # Calculate the final result
    expected_max_magnitude = pi_val / (two_val * log2_val)
    
    # Print the final equation and the values of its components as requested
    print("The final equation for the expected maximum magnitude is E[M] = pi / (2 * log(2))")
    print(f"pi = {pi_val}")
    print(f"2 = {two_val}")
    print(f"log(2) = {log2_val}")
    
    # Print the calculation steps with the numerical values
    print("\nCalculating the result:")
    print(f"E[M] = {pi_val} / ({two_val} * {log2_val})")
    
    # Print the final answer
    print(f"The expected maximum earthquake magnitude is: {expected_max_magnitude}")

# Execute the function
solve_earthquake_magnitude()