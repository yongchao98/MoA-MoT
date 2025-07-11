def display_potential_formula():
    """
    This function prints the derived formula for the electric potential
    Phi(x, y) in the region 0 <= y <= a, breaking it down into its
    numerator and denominator for clarity.
    """

    # Define strings to represent the symbolic components of the formula.
    sigma_0_str = "sigma_0"
    k_str = "k"
    x_str = "x"
    y_str = "y"
    a_str = "a"
    b_str = "b"
    epsilon_1_str = "epsilon_1"
    epsilon_2_str = "epsilon_2"

    # Construct the string for the numerator of the formula.
    numerator = f"-{sigma_0_str} * sinh({k_str}*{b_str}) * sinh({k_str}*({y_str} - {a_str})) * sin({k_str}*{x_str})"

    # Construct the string for the denominator of the formula.
    denominator_part1 = f"{epsilon_2_str}*cosh({k_str}*{a_str})*sinh({k_str}*{b_str})"
    denominator_part2 = f"{epsilon_1_str}*sinh({k_str}*{a_str})*cosh({k_str}*{b_str})"
    denominator = f"{k_str} * [{denominator_part1} + {denominator_part2}]"

    print("The derived electric potential Phi(x, y) in the region 0 <= y <= a is:")
    print("-" * 70)
    print(f"Phi(x, y) = ({numerator}) / ({denominator})")
    print("-" * 70)
    print("\nThis formula corresponds to the expression for the region 0 < y < a in answer choice A.")

if __name__ == '__main__':
    display_potential_formula()