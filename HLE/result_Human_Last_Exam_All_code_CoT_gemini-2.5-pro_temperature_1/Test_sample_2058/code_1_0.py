from fractions import Fraction

def solve():
    """
    This program calculates the total mass M(A_5, rho, 2).
    The formula for the total mass is:
    M(G, rho, p) = (1/|G|) * sum_{K} |Aut(Gal(K/Q_p))| * |Inj(Gal(K/Q_p), G)| * p^(-c_rho(K))
    where the sum is over isomorphism classes of Galois extensions K/Q_p.

    G = A_5, |G| = 60
    p = 2
    rho is the 5D permutation representation. The conductor c_rho(K) = c_{chi_4}(K),
    where chi_4 is the 4D irreducible character of A_5.

    We sum the contributions from each possible solvable Galois group H <= A_5.
    Contribution(H) = (1/|G|) * |Aut(H)| * |Inj(H, A_5)| * sum_{K, Gal=H} 2^(-c_rho(K))
    """

    total_mass = Fraction(0)
    G_order = 60

    # H = {1}, the trivial group
    # One extension K=Q_2, c_rho(K)=0.
    # |Aut({1})| = 1, |Inj({1}, A_5)| = 1
    # Sum_K = 2^0 = 1
    H1_aut = 1
    H1_inj = 1
    H1_sum_K = 1
    mass_H1 = Fraction(H1_aut * H1_inj * H1_sum_K, G_order)
    total_mass += mass_H1
    print(f"Contribution from H = C1: {mass_H1}")

    # H = C_2 (cyclic group of order 2)
    # There are 7 quadratic extensions of Q_2.
    # c_rho(K) = 2 * c(K).
    # 1 unramified ext (c(K)=0), c_rho = 0.
    # 2 ramified exts (c(K)=2), c_rho = 4.
    # 4 ramified exts (c(K)=3), c_rho = 6.
    # |Aut(C_2)| = 1, |Inj(C_2, A_5)| = 15.
    H2_aut = 1
    H2_inj = 15
    H2_sum_K = (Fraction(1, 2**0) * 1 +
                Fraction(1, 2**4) * 2 +
                Fraction(1, 2**6) * 4)
    mass_H2 = Fraction(H2_aut * H2_inj, G_order) * H2_sum_K
    total_mass += mass_H2
    print(f"Contribution from H = C2: {mass_H2}")

    # H = C_3
    # There are 2 cubic extensions of Q_2 with group C_3.
    # c_rho(K) = 2 * c(K).
    # 1 unramified ext (c(K)=0), c_rho = 0.
    # 1 tamely ramified ext (c(K)=1), c_rho = 2.
    # |Aut(C_3)| = 2, |Inj(C_3, A_5)| = 20.
    H3_aut = 2
    H3_inj = 20
    H3_sum_K = (Fraction(1, 2**0) * 1 +
                Fraction(1, 2**2) * 1)
    mass_H3 = Fraction(H3_aut * H3_inj, G_order) * H3_sum_K
    total_mass += mass_H3
    print(f"Contribution from H = C3: {mass_H3}")

    # H = V_4 = C_2 x C_2
    # There are 7 V_4-extensions of Q_2.
    # c_rho(K) = c(K1)+c(K2)+c(K3) for the three quadratic subextensions.
    # Conductor sums: one is 4, two are 6, four are 8.
    # |Aut(V_4)| = GL(2,F_2) = 6. |Inj(V_4, A_5)| = 30.
    H4_aut = 6
    H4_inj = 30
    H4_sum_K = (Fraction(1, 2**4) * 1 +
                Fraction(1, 2**6) * 2 +
                Fraction(1, 2**8) * 4)
    mass_H4 = Fraction(H4_aut * H4_inj, G_order) * H4_sum_K
    total_mass += mass_H4
    print(f"Contribution from H = V4: {mass_H4}")

    # H = S_3
    # There are 2 S_3-extensions of Q_2 with inertia group A_3.
    # For these, c_rho = 2.
    # |Aut(S_3)| = 6. |Inj(S_3, A_5)| = 60.
    H5_aut = 6
    H5_inj = 60
    H5_sum_K = Fraction(1, 2**2) * 2
    mass_H5 = Fraction(H5_aut * H5_inj, G_order) * H5_sum_K
    total_mass += mass_H5
    print(f"Contribution from H = S3: {mass_H5}")

    # H = D_10 (dihedral of order 10)
    # There is 1 D_10-extension, which is tamely ramified with inertia C_5.
    # c_rho = 4.
    # |Aut(D_10)| = 20. |Inj(D_10, A_5)| = 120.
    H6_aut = 20
    H6_inj = 120
    H6_sum_K = Fraction(1, 2**4) * 1
    mass_H6 = Fraction(H6_aut * H6_inj, G_order) * H6_sum_K
    total_mass += mass_H6
    print(f"Contribution from H = D10: {mass_H6}")

    # H = A_4
    # There are 2 A_4-extensions of Q_2.
    # For these, c_rho = 3 for one and c_rho = 4 for the other.
    # |Aut(A_4)| = 24. |Inj(A_4, A_5)| = 120.
    H7_aut = 24
    H7_inj = 120
    H7_sum_K = (Fraction(1, 2**3) * 1 +
                Fraction(1, 2**4) * 1)
    mass_H7 = Fraction(H7_aut * H7_inj, G_order) * H7_sum_K
    total_mass += mass_H7
    print(f"Contribution from H = A4: {mass_H7}")

    print(f"\nTotal Mass Calculation:")
    c1_term = mass_H1
    c2_term = mass_H2
    c3_term = mass_H3
    v4_term = mass_H4
    s3_term = mass_H5
    d10_term = mass_H6
    a4_term = mass_H7
    
    final_sum = c1_term + c2_term + c3_term + v4_term + s3_term + d10_term + a4_term
    
    print(f"M = {c1_term} (C1) + {c2_term} (C2) + {c3_term} (C3) + {v4_term} (V4) + {s3_term} (S3) + {d10_term} (D10) + {a4_term} (A4)")
    print(f"M = {c1_term.numerator}/{c1_term.denominator} + {c2_term.numerator}/{c2_term.denominator} + {c3_term.numerator}/{c3_term.denominator} + {v4_term.numerator}/{v4_term.denominator} + {s3_term.numerator}/{s3_term.denominator} + {d10_term.numerator}/{d10_term.denominator} + {a4_term.numerator}/{a4_term.denominator}")
    
    # Simplify terms for display
    combined_c2_v4 = c2_term + v4_term
    combined_c1_c3 = c1_term + c3_term
    combined_s3_d10_a4 = s3_term + d10_term + a4_term
    
    print(f"Combining terms with common denominators:")
    print(f"C2 + V4 terms = {c2_term} + {v4_term} = {combined_c2_v4} = {combined_c2_v4.numerator}/{combined_c2_v4.denominator}")
    print(f"C1 + C3 terms = {c1_term} + {c3_term} = {combined_c1_c3} = {combined_c1_c3.numerator}/{combined_c1_c3.denominator}")
    print(f"S3 + D10 + A4 terms = {s3_term} + {d10_term} + {a4_term} = {combined_s3_d10_a4} = {combined_s3_d10_a4.numerator}/{combined_s3_d10_a4.denominator}")

    # Summing up simplified terms
    final_sum_simplified = combined_c1_c3 + combined_c2_v4 + combined_s3_d10_a4
    print(f"Summing these parts: {combined_c1_c3} + {combined_c2_v4} + {combined_s3_d10_a4} = {final_sum_simplified}")
    
    print(f"\nThe total mass is {total_mass.numerator}/{total_mass.denominator}")

solve()