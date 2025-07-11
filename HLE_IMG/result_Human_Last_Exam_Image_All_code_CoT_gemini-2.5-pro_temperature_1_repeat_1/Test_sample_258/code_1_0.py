import math

def calculate_cp_length():
    """
    This function calculates the length of segment CP given the length of AB.
    The derivation shows that CP = AB * φ, where φ is the golden ratio.
    """
    # The given length of the segment AB
    ab_length = 18468

    # The formula for the golden ratio is (1 + sqrt(5)) / 2
    sqrt_5 = math.sqrt(5)
    golden_ratio = (1 + sqrt_5) / 2

    # Calculate the length of CP
    cp_length = ab_length * golden_ratio

    # Print the equation and the steps of the calculation
    print(f"Given the length of AB = {ab_length}")
    print("The derived formula for the length of CP is: CP = AB * (1 + sqrt(5)) / 2")
    print("\n--- Calculation Steps ---")
    print(f"CP = {ab_length} * (1 + {sqrt_5}) / 2")
    print(f"CP = {ab_length} * {golden_ratio}")
    print(f"CP = {cp_length}")

    # Print the final answer rounded to 4 decimal points
    print(f"\nThe final length of CP rounded to 4 decimal places is: {cp_length:.4f}")

calculate_cp_length()