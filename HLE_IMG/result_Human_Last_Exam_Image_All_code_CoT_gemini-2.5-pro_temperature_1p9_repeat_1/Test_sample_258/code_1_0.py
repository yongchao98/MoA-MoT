import math

def solve_for_cp():
    """
    Calculates the length of CP given the length of AB, based on the
    geometric construction provided.
    """
    # Given length of segment AB
    ab_length = 18468

    # The derived formula is CP = AB * (sqrt(5) + 1) / 2
    # This is AB times the golden ratio.
    sqrt_5 = math.sqrt(5)
    golden_ratio = (sqrt_5 + 1) / 2
    
    cp_length = ab_length * golden_ratio

    # Outputting the final equation with the numbers plugged in, as requested.
    print(f"Given AB = {ab_length}")
    print("The formula is: CP = AB * (sqrt(5) + 1) / 2")
    print("\nCalculation:")
    # Printing each number in the final equation.
    print(f"CP = {ab_length} * ({sqrt_5} + 1) / 2")
    print(f"CP = {ab_length} * {sqrt_5 + 1} / 2")
    print(f"CP = {ab_length * (sqrt_5 + 1)} / 2")
    print(f"CP = {cp_length}")

    # Final answer formatted to 4 decimal places
    print(f"\nThe length of CP is {cp_length:.4f}")

solve_for_cp()