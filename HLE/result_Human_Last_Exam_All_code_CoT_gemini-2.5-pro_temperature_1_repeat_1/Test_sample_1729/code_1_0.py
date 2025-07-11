def solve_and_print_formula():
    """
    This function calculates and prints the symbolic formula for the probability P_m.
    """

    # Numerator and denominator components of the probability P_m
    # These are the integer constants that define the formula.
    numerator_m_1 = 3
    
    numerator_m_ge_2_coeff_m = 3
    numerator_m_ge_2_const = 1
    
    denominator_part1_coeff_m = 2
    denominator_part1_const = 1
    
    denominator_part2_coeff_m = 4
    denominator_part2_const = 1

    print("The probability P_m is given by a piecewise formula depending on the value of m.")
    print("-" * 60)

    # Case m = 1
    print("Case 1: m = 1")
    denom_1_val = (denominator_part1_coeff_m * 1 + denominator_part1_const) * (denominator_part2_coeff_m * 1 + denominator_part2_const)
    print(f"The number of favorable pairs (i,j) is {numerator_m_1}.")
    print(f"The total number of pairs is ({denominator_part1_coeff_m}*1 + {denominator_part1_const}) * ({denominator_part2_coeff_m}*1 + {denominator_part2_const}) = {denom_1_val}.")
    print(f"The equation for P_1 is:")
    print(f"    P_1 = {numerator_m_1} / {denom_1_val} = 1/5")
    print("-" * 60)

    # Case m >= 2
    print("Case 2: m >= 2")
    print(f"The number of favorable pairs (i,j) is {numerator_m_ge_2_coeff_m}*m + {numerator_m_ge_2_const}.")
    print(f"The total number of pairs is ({denominator_part1_coeff_m}*m + {denominator_part1_const})*({denominator_part2_coeff_m}*m + {denominator_part2_const}).")
    print(f"The equation for P_m is:")
    print(f"    P_m = ({numerator_m_ge_2_coeff_m}*m + {numerator_m_ge_2_const}) / (({denominator_part1_coeff_m}*m + {denominator_part1_const})*({denominator_part2_coeff_m}*m + {denominator_part2_const}))")
    print("-" * 60)

# Execute the function to display the result
solve_and_print_formula()