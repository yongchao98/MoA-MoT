def find_pm_formula():
    """
    This function states the final formula for the probability P_m,
    based on the parity of the positive integer m.
    The derivation finds that the number of successful (i,j) pairs depends
    on whether m is even or odd, while the total number of pairs is always
    the same for a given m.
    """
    
    # The total number of ways to choose i and j from {1, ..., 4m+2} with i < j is C(4m+2, 2).
    # C(4m+2, 2) = (4m+2)(4m+1)/2 = (2m+1)(4m+1).
    
    # The number of successful pairs, K_m, is found to be 3 if m is odd, and 4 if m is even.
    
    print("The formula for P_m depends on the parity of m.")
    
    print("\nIf m is an odd positive integer:")
    numerator_odd = 3
    denominator_part1_coeff = 2
    denominator_part1_const = 1
    denominator_part2_coeff = 4
    denominator_part2_const = 1
    
    print("The number of successful (i,j) pairs is 3.")
    print("The total number of pairs is (2*m + 1)*(4*m + 1).")
    print("The probability P_m is given by the formula:")
    print(f"P_m = {numerator_odd} / (({denominator_part1_coeff}*m + {denominator_part1_const}) * ({denominator_part2_coeff}*m + {denominator_part2_const}))")

    print("\nIf m is an even positive integer:")
    numerator_even = 4
    # Denominator parts are the same
    
    print("The number of successful (i,j) pairs is 4.")
    print("The total number of pairs is (2*m + 1)*(4*m + 1).")
    print("The probability P_m is given by the formula:")
    print(f"P_m = {numerator_even} / (({denominator_part1_coeff}*m + {denominator_part1_const}) * ({denominator_part2_coeff}*m + {denominator_part2_const}))")

find_pm_formula()