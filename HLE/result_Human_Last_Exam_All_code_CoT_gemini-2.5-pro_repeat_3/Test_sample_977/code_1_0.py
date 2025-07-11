def print_potential_formula():
    """
    This function prints the derived formula for the electric potential Phi(x, y).
    The derivation is based on solving Laplace's equation with the given boundary conditions.
    """
    # Define symbolic variables as strings for printing the formula
    sigma_0 = "sigma_0"
    k = "k"
    x = "x"
    y = "y"
    a = "a"
    b = "b"
    epsilon_1 = "epsilon_1"
    epsilon_2 = "epsilon_2"

    # Common denominator term derived from boundary conditions
    denominator = f"k * [ {epsilon_2}*cosh({k}*{a})*sinh({k}*{b}) + {epsilon_1}*sinh({k}*{a})*cosh({k}*{b}) ]"

    # Numerator for the region 0 < y < a
    numerator_2 = f"-{sigma_0} * sinh({k}*{b}) * sinh({k}*({y} - {a})) * sin({k}*{x})"
    
    # Numerator for the region -b < y < 0
    numerator_1 = f"{sigma_0} * sinh({k}*{a}) * sinh({k}*({y} + {b})) * sin({k}*{x})"

    print("The electric potential Phi(x, y) is given by a piecewise function corresponding to Answer Choice A:")
    print("="*80)

    # Print the potential for the region 0 < y < a
    print("For the region 0 < y < a:")
    print("Phi(x, y) =")
    print(f"  {numerator_2}")
    print("  " + "-" * (len(numerator_2)-2))
    print(f"  {denominator}\n")

    # Print the potential for the region -b < y < 0
    print("For the region -b < y < 0:")
    print("Phi(x, y) =")
    print(f"  {numerator_1}")
    print("  " + "-" * (len(numerator_1)-2))
    print(f"  {denominator}")
    
    print("="*80)
    print("The question asks for the potential in the region 0 <= y <= a, which is the first expression above.")

if __name__ == '__main__':
    print_potential_formula()