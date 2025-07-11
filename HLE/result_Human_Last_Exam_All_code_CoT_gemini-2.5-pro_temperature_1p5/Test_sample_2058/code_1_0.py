def solve():
    """
    This function calculates the total mass M(A_5, rho, 2) based on the detailed breakdown.
    The total mass is the sum of contributions from each subgroup type.
    """
    
    # Contribution from C_1
    c1_contr = 1/60
    
    # Contribution from C_2
    c2_contr = 5/8
    
    # Contribution from C_3
    c3_contr = 1/3
    
    # Contribution from C_5
    c5_contr = 2/5
    
    # Contribution from V_4
    c_v4 = 105/8192
    
    # Contribution from S_3
    c_s3 = 3/8
    
    # Contribution from D_10
    c_d10 = 1/8
    
    # Contribution from A_4
    c_a4 = 3/8
    
    # Sum of rational numbers without the 8192 denominator
    from fractions import Fraction
    total_mass_simple_fracs = Fraction(c1_contr) + Fraction(c2_contr) + Fraction(c3_contr) + Fraction(c5_contr) + Fraction(c_s3) + Fraction(c_d10) + Fraction(c_a4)
    
    total_mass = total_mass_simple_fracs + Fraction(c_v4)
    
    numerator = total_mass.numerator
    denominator = total_mass.denominator

    print(f"The total mass M(A_5, rho, 2) is the sum of contributions from all solvable subgroups.")
    print(f"Contribution from C1: 1/60")
    print(f"Contribution from C2: 5/8")
    print(f"Contribution from C3: 1/3")
    print(f"Contribution from C5: 2/5")
    print(f"Contribution from V4: 105/8192")
    print(f"Contribution from S3: 3/8")
    print(f"Contribution from D10: 1/8")
    print(f"Contribution from A4: 3/8")
    print(f"Summing these values: M = 1/60 + 5/8 + 1/3 + 2/5 + 105/8192 + 3/8 + 1/8 + 3/8")
    # Simplify part of the sum
    sum_of_eighths = Fraction(5,8) + Fraction(3,8) + Fraction(1,8) + Fraction(3,8)
    sum_others = Fraction(1,60) + Fraction(1,3) + Fraction(2,5)
    
    # print intermediate sum steps
    print(f"Sum of eighths: 5/8 + 3/8 + 1/8 + 3/8 = {sum_of_eighths.numerator}/{sum_of_eighths.denominator} = 3/2")
    print(f"Sum of others: 1/60 + 1/3 + 2/5 = {sum_others.numerator}/{sum_others.denominator} = 3/4")
    print(f"Intermediate sum: 3/4 + 3/2 = 9/4")
    print(f"Final sum: M = 9/4 + 105/8192")
    print(f"In common denominator: M = (9 * 2048) / 8192 + 105/8192 = 18432/8192 + 105/8192")
    
    final_numerator = 18432 + 105
    final_denominator = 8192
    
    print(f"Final Answer: {final_numerator}/{final_denominator}")

solve()