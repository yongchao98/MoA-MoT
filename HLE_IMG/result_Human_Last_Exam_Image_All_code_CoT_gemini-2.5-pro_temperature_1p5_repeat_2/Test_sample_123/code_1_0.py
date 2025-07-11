from fractions import Fraction

def solve_bran_castle_challenge():
    """
    This function sets up and solves the system of linear equations for the
    Bran Castle escape challenge to find the probability of reaching the Treasure Room.
    """
    # We solve the 2x2 system for p_gh and p_kh derived from the full system.
    # The system is:
    # 449 * p_gh - 308 * p_kh = 69
    # -11 * p_gh +  14 * p_kh = 3

    a1, b1, c1 = Fraction(449), Fraction(-308), Fraction(69)
    a2, b2, c2 = Fraction(-11), Fraction(14), Fraction(3)

    # Using Cramer's rule to solve the 2x2 system
    D = a1 * b2 - b1 * a2
    D_gh = c1 * b2 - b1 * c2
    D_kh = a1 * c2 - c1 * a2

    p_gh = D_gh / D
    p_kh = D_kh / D

    # Now calculate p_me using the relation: p_me = (2 * p_gh + p_kh) / 3
    # This relation is derived from:
    # (1) p_me = 1/2 * p_gh + 1/2 * p_ck
    # (3) p_ck = 1/2 * p_me + 1/2 * p_kh
    # Substituting (3) into (1) gives the relation.

    p_me_numerator = 2 * p_gh + p_kh
    p_me = p_me_numerator / 3

    print("Step 1: Solve for the key probabilities p_gh (Great Hall) and p_kh (Knights' Hall).")
    print(f"Probability from Great Hall (p_gh) = {p_gh.numerator}/{p_gh.denominator}")
    print(f"Probability from Knights' Hall (p_kh) = {p_kh.numerator}/{p_kh.denominator}")
    print("\nStep 2: Calculate the final probability for the Main Entrance (p_me).")
    print("The formula is: p_me = (2 * p_gh + p_kh) / 3\n")
    print("Substituting the values:")
    print(f"p_me = (2 * ({p_gh.numerator}/{p_gh.denominator}) + {p_kh.numerator}/{p_kh.denominator}) / 3")
    
    # Show the intermediate calculation with a common denominator
    term1 = 2 * p_gh
    common_denominator = p_kh.denominator
    term1_num_new = term1.numerator * (common_denominator // term1.denominator)

    print(f"p_me = (({term1_num_new}/{common_denominator}) + ({p_kh.numerator}/{p_kh.denominator})) / 3")
    
    p_me_num_val = p_me_numerator.numerator
    p_me_den_val = p_me_numerator.denominator
    print(f"p_me = (({p_me_num_val}/{p_me_den_val})) / 3")

    final_num = p_me.numerator
    final_den = p_me.denominator
    print(f"p_me = {p_me_num_val} / ({p_me_den_val} * 3)")
    print(f"\nThe final probability as an irreducible fraction is: {final_num}/{final_den}")

solve_bran_castle_challenge()
print("\n<<<109/161>>>")