from fractions import Fraction

def calculate_total_mass():
    """
    This script calculates the total mass M(A_5, rho, 2) by summing the contributions
    from each relevant subgroup of A_5.
    """
    
    # Contribution from H = {e}
    # N(H) = 60. One extension K with c=0.
    # Term = 1/60 * 2^0
    term_e = Fraction(1, 60)
    
    # Contribution from H = C_2
    # N(H) = 4. There are 7 extensions.
    # 1 unramified (c=0), 2 with conductor c=4, 4 with conductor c=6.
    # Term = 1/4 * (1*2^-0 + 2*2^-4 + 4*2^-6)
    sum_C2 = Fraction(1, 1) + 2 * Fraction(1, 16) + 4 * Fraction(1, 64)
    term_C2 = Fraction(1, 4) * sum_C2
    
    # Contribution from H = C_3
    # N(H) = 6. There are 2 extensions.
    # 1 unramified (c=0), 1 tamely ramified (c=2).
    # Term = 1/6 * (2^-0 + 2^-2)
    sum_C3 = Fraction(1, 1) + Fraction(1, 4)
    term_C3 = Fraction(1, 6) * sum_C3
    
    # Contribution from H = V_4
    # N(H) = 12. There are 19 extensions.
    # The sum of 2^(-c) over these extensions is 139/1024.
    # Term = 1/12 * (139/1024)
    sum_V4 = Fraction(1, 16) + 2*Fraction(1,64) + 8*Fraction(1,256) + 4*Fraction(1,512) + 2*Fraction(1,1024) + 2*Fraction(1,2048)
    # sum_V4 simplifies to Fraction(139, 1024)
    term_V4 = Fraction(1, 12) * sum_V4
    
    # Contribution from H = A_4
    # N(H) = 12. There is 1 extension.
    # Conductor c=4.
    # Term = 1/12 * 2^-4
    term_A4 = Fraction(1, 12) * Fraction(1, 16)

    # Total Mass is the sum of all contributions.
    total_mass = term_e + term_C2 + term_C3 + term_V4 + term_A4

    # Print the equation
    print("The total mass M is the sum of the following contributions:")
    print(f"M = {term_e.numerator}/{term_e.denominator} ({'H={e}'})")
    print(f"  + {term_C2.numerator}/{term_C2.denominator} ({'H=C2'})")
    print(f"  + {term_C3.numerator}/{term_C3.denominator} ({'H=C3'})")
    print(f"  + {term_V4.numerator}/{term_V4.denominator} ({'H=V4'})")
    print(f"  + {term_A4.numerator}/{term_A4.denominator} ({'H=A4'})")

    print(f"\nThe equation with common denominator {total_mass.denominator} is:")
    
    # To show the equation with a common denominator
    common_denom = total_mass.denominator
    
    t1 = term_e * common_denom / term_e.denominator
    t2 = term_C2 * common_denom / term_C2.denominator
    t3 = term_C3 * common_denom / term_C3.denominator
    t4 = term_V4 * common_denom / term_V4.denominator
    t5 = term_A4 * common_denom / term_A4.denominator

    print(f"M = ({int(t1)} + {int(t2)} + {int(t3)} + {int(t4)} + {int(t5)}) / {common_denom}")
    print(f"M = {int(t1)+int(t2)+int(t3)+int(t4)+int(t5)} / {common_denom}")
    
    print(f"\nThe final result as a fraction in lowest terms:")
    print(f"M = {total_mass.numerator}/{total_mass.denominator}")

calculate_total_mass()