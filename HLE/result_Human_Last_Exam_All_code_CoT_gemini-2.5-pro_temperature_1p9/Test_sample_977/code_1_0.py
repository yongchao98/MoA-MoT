def solve_potential_problem():
    """
    This function formulates and prints the final expression for the electric potential
    in the region 0 < y < a based on the step-by-step derivation.
    """
    
    # Define symbolic variables as strings for presentation
    sigma_0 = "sigma_0"
    k = "k"
    x = "x"
    y = "y"
    a = "a"
    b = "b"
    epsilon_1 = "epsilon_1"
    epsilon_2 = "epsilon_2"

    # Construct the numerator of the expression
    # Numerator = -sigma_0 * sinh(kb) * sinh(k(y-a)) * sin(kx)
    numerator = f"-{sigma_0} * sinh({k}{b}) * sinh({k}({y} - {a})) * sin({k}{x})"
    
    # Construct the denominator of the expression
    # Denominator = k * [epsilon_2 * cosh(ka) * sinh(kb) + epsilon_1 * sinh(ka) * cosh(kb)]
    denominator_part1 = f"{epsilon_2}*cosh({k}{a})*sinh({k}{b})"
    denominator_part2 = f"{epsilon_1}*sinh({k}{a})*cosh({k}{b})"
    denominator = f"{k} * [ {denominator_part1} + {denominator_part2} ]"

    # Print the final result in a readable format
    print("The electric potential Phi(x, y) in the region 0 < y < a is:")
    print("\n")
    print(f"Phi(x, y) = ({numerator}) / ({denominator})")
    print("\nThis matches option A.")


if __name__ == "__main__":
    solve_potential_problem()
