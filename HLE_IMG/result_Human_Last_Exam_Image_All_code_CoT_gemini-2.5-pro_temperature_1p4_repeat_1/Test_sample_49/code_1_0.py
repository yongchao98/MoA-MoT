import math

def solve_cutoff_frequency():
    """
    This function derives and prints the formula for the cutoff frequency
    of the given infinite ladder network.
    """
    
    # The final derived equation for the cutoff frequency is:
    # omega_c = 1 / ((1 + sqrt(3)) * r * C)
    
    # Extract the numerical values from the formula to print them as requested.
    numerator = 1
    term1_in_denominator = 1
    term2_in_denominator_sqrt = 3
    
    print("The derived Thévenin equivalent resistance is: R_th = (1 + sqrt(3)) * r")
    print("\nThe formula for the angular cutoff frequency is: omega_c = 1 / (R_th * C)")
    print("\nSubstituting R_th, the final equation is:")
    
    # Print the final equation with each number explicitly shown.
    print(f"omega_c = {numerator} / (({term1_in_denominator} + sqrt({term2_in_denominator_sqrt})) * r * C)")

solve_cutoff_frequency()