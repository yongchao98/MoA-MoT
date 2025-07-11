import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of segment CP based on the length of AB.

    The derivation shows that the ratio of CP to AB is the golden ratio, phi.
    CP = AB * phi
    phi = (1 + sqrt(5)) / 2
    """
    
    # Calculate the value of sqrt(5)
    sqrt_5 = math.sqrt(5)
    
    # Calculate the golden ratio, phi
    phi = (1 + sqrt_5) / 2
    
    # Calculate the length of CP
    cp_length = ab_length * phi
    
    # Print the explanation and the equation with values
    print(f"Given the length of segment AB = {ab_length}")
    print("The length of segment CP can be found using the golden ratio (phi).")
    print("The formula is: CP = AB * (1 + sqrt(5)) / 2")
    print("\nSubstituting the values into the equation:")
    print(f"CP = {ab_length} * (1 + {sqrt_5:.4f}) / 2")
    print(f"CP = {ab_length} * ({1 + sqrt_5:.4f}) / 2")
    print(f"CP = {ab_length} * {phi:.4f}")
    
    # Print the final result formatted to 4 decimal places
    print(f"\nThe final calculated length of CP is: {cp_length:.4f}")
    
    return cp_length

# The given length of AB
ab_length_input = 18468
result = calculate_cp_length(ab_length_input)
print(f"\n<<<{result:.4f}>>>")
