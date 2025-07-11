import math

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet based on the simplified model.
    """
    # Given constants
    sigma_0 = 7.43e-7  # in e/nm
    R_0 = 30.0         # in nm

    # The formula for the total charge simplifies to Q = 4 * pi * sigma_0 * R_0
    # The complex terms and perturbation parameters cancel out or integrate to zero.
    
    # Calculate the total charge
    Q = 4 * math.pi * sigma_0 * R_0

    # Print the steps of the calculation as requested
    print("Based on the simplification, the total charge Q is calculated using the formula:")
    print("Q = 4 * pi * sigma_0 * R_0")
    print("\nSubstituting the given values:")
    print(f"sigma_0 = {sigma_0} e/nm")
    print(f"R_0 = {R_0} nm")
    print(f"pi = {math.pi}")
    print("\nThe calculation is:")
    print(f"Q = 4 * {math.pi} * {sigma_0} * {R_0}")
    
    # Print the final result
    print(f"\nThe total charge Q is: {Q} e")

calculate_total_charge()