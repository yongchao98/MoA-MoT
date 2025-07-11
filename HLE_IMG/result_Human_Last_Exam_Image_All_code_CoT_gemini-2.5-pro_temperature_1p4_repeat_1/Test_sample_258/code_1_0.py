import math

def solve_geometry():
    """
    Calculates the length of CP given the length of AB based on the described
    geometric construction.
    """
    # The given length of the segment AB.
    ab_length = 18468

    # The geometric analysis shows that the ratio of the length of CP to the length of AB
    # is the golden ratio, phi (φ).
    # Formula: CP = AB * (1 + sqrt(5)) / 2
    
    # Value of sqrt(5)
    sqrt_5_val = math.sqrt(5)
    
    # The golden ratio, phi
    phi = (1 + sqrt_5_val) / 2
    
    # Calculate the length of CP
    cp_length = ab_length * phi
    
    # Print the steps of the final calculation as requested
    print("The relationship found is CP = AB * φ, where φ is the golden ratio.")
    print("The formula for the golden ratio is: (1 + sqrt(5)) / 2")
    print("\nCalculating the length of CP:")
    print(f"CP = {ab_length} * (1 + {sqrt_5_val}) / 2")
    print(f"CP = {ab_length} * ({1 + sqrt_5_val}) / 2")
    print(f"CP = {ab_length * (1 + sqrt_5_val)} / 2")
    
    # Print the final answer, formatted to 4 decimal places.
    print(f"\nThe final length of CP is: {cp_length:.4f}")

solve_geometry()