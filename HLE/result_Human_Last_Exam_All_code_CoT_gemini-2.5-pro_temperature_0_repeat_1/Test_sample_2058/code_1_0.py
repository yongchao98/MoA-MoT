from fractions import Fraction

def calculate_mass():
    """
    Calculates the total mass M(A_5, rho, 2).
    """
    # Order of the group A_5
    order_A5 = 60

    # Contributions from each type of solvable subgroup H
    # Term = N_H * S_H where N_H is the number of subgroups and
    # S_H is the sum over extensions with Galois group H.

    # H = C_1 (trivial group)
    N_C1 = 1
    S_C1 = Fraction(1)
    term_C1 = N_C1 * S_C1

    # H = C_2
    N_C2 = 15
    # 1 unramified ext (f=0), 2 ramified (f=2), 4 ramified (f=3)
    # c_rho = 3*f(L)
    S_C2 = Fraction(1, 2**(3*0)) + 2 * Fraction(1, 2**(3*2)) + 4 * Fraction(1, 2**(3*3))
    term_C2 = N_C2 * S_C2

    # H = C_3
    N_C3 = 10
    # 1 unramified ext (c_rho=0), 1 ramified (c_rho=3)
    S_C3 = Fraction(1, 2**0) + Fraction(1, 2**3)
    term_C3 = N_C3 * S_C3

    # H = V_4
    N_V4 = 5
    # 2 extensions with c_rho = 6 and 9
    S_V4 = Fraction(1, 2**6) + Fraction(1, 2**9)
    term_V4 = N_V4 * S_V4

    # H = S_3
    N_S3 = 10
    # 5 extensions with c_rho = 5, 7, 9, 12, 13
    S_S3 = (Fraction(1, 2**5) + Fraction(1, 2**7) + Fraction(1, 2**9) +
            Fraction(1, 2**12) + Fraction(1, 2**13))
    term_S3 = N_S3 * S_S3

    # Total sum inside the formula
    total_sum = term_C1 + term_C2 + term_C3 + term_V4 + term_S3

    # Final mass
    mass = total_sum / order_A5

    # Print the equation
    print(f"M(A_5, rho, 2) = 1/60 * ( (1 * {S_C1}) + (15 * {S_C2}) + (10 * {S_C3}) + (5 * {S_V4}) + (10 * {S_S3}) )")
    print(f"M(A_5, rho, 2) = 1/60 * ( {term_C1} + {term_C2} + {term_C3} + {term_V4} + {term_S3} )")
    print(f"M(A_5, rho, 2) = 1/60 * ( {total_sum} )")
    print(f"M(A_5, rho, 2) = {mass}")

calculate_mass()