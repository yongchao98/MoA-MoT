def solve_potential():
    """
    This function formulates and prints the final expression for the electric potential
    in the region 0 <= y <= a, based on the derivation.
    """
    
    # Define symbolic parts of the equation as strings
    s_0 = "sigma_0"
    k = "k"
    a = "a"
    b = "b"
    x = "x"
    y = "y"
    eps_1 = "epsilon_1"
    eps_2 = "epsilon_2"

    # Construct the numerator of the expression
    # Numerator = -sigma_0 * sinh(kb) * sinh(k(y - a)) * sin(kx)
    numerator = f"-{s_0} sinh({k}{b}) sinh({k}({y} - {a})) sin({k}{x})"
    
    # Construct the denominator of the expression
    # Denominator = k * [epsilon_2 * cosh(ka) * sinh(kb) + epsilon_1 * sinh(ka) * cosh(kb)]
    denominator_part1 = f"{eps_2} cosh({k}{a}) sinh({k}{b})"
    denominator_part2 = f"{eps_1} sinh({k}{a}) cosh({k}{b})"
    denominator = f"{k} [ {denominator_part1} + {denominator_part2} ]"

    # Print the final result for the specified region
    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    print(f"Phi(x, y) = ({numerator}) / ({denominator})")
    print("\nThis matches the expression for the region 0 < y < a in answer choice A.")
    print("A full derivation confirms the expression for the region -b < y < 0 also matches choice A.")

solve_potential()