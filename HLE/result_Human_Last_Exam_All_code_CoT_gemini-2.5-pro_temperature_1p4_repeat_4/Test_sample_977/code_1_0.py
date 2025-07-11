import math

def print_potential_formula():
    """
    Prints the symbolic formula for the electric potential Phi(x, y).
    The formula is broken down into its constituent parts for clarity.
    """

    print("The derived electric potential Phi(x, y) has two expressions for the two regions:")

    # Common Denominator
    denominator = "k * [epsilon_2 * cosh(ka) * sinh(kb) + epsilon_1 * sinh(ka) * cosh(kb)]"
    
    print("\n--- For the region 0 < y < a ---")
    
    # Numerator for Region 2 (0 < y < a)
    num_2_part_1 = "-sigma_0"
    num_2_part_2 = "sinh(kb)"
    num_2_part_3 = "sinh(k*(y - a))"
    num_2_part_4 = "sin(kx)"
    numerator_2 = f"{num_2_part_1} * {num_2_part_2} * {num_2_part_3} * {num_2_part_4}"
    
    print("Potential Phi_2(x, y):")
    print(f"  Numerator:   {numerator_2}")
    print(f"  Denominator: {denominator}")
    
    print("\n--- For the region -b < y < 0 ---")
    
    # Numerator for Region 1 (-b < y < 0)
    num_1_part_1 = "sigma_0"
    num_1_part_2 = "sinh(ka)"
    num_1_part_3 = "sinh(k*(y + b))"
    num_1_part_4 = "sin(kx)"
    numerator_1 = f"{num_1_part_1} * {num_1_part_2} * {num_1_part_3} * {num_1_part_4}"

    print("Potential Phi_1(x, y):")
    print(f"  Numerator:   {numerator_1}")
    print(f"  Denominator: {denominator}")
    
    print("\nThis combined expression for Phi(x, y) corresponds to Answer Choice A.")

# Execute the function to print the formula
print_potential_formula()