import math

def find_cp_length(ab_length):
    """
    Calculates the length of CP given the length of AB based on the geometric construction.
    The derived relationship is CP = AB * φ, where φ is the golden ratio.
    """
    
    # The golden ratio, φ = (1 + sqrt(5)) / 2
    golden_ratio = (1 + math.sqrt(5)) / 2
    
    # Calculate the length of CP
    cp_length = ab_length * golden_ratio
    
    # Print the equation with the given numbers for clarity
    print(f"Given AB = {ab_length}")
    print("The derived formula is: CP = AB * (1 + √5) / 2")
    print(f"Plugging in the numbers: CP = {ab_length} * (1 + {math.sqrt(5):.9f}) / 2")
    print(f"CP = {ab_length} * {golden_ratio:.9f}")
    
    # Print the final result formatted to 4 decimal places
    print(f"\nThe calculated length of CP is: {cp_length:.4f}")

# The given length of the segment AB
ab_value = 18468
find_cp_length(ab_value)