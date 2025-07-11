import itertools

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    res = []
    for i in range(max(len(p1), len(p2))):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res.append((c1 + c2) % 2)
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return tuple(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if not p1 or not p2 or p1 == (0,) or p2 == (0,):
        return (0,)
    res = [0] * (len(p1) + len(p2) - 1)
    for i1, c1 in enumerate(p1):
        for i2, c2 in enumerate(p2):
            res[i1 + i2] = (res[i1 + i2] + c1 * c2) % 2
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return tuple(res)

def poly_deg(p):
    """Calculates the degree of a polynomial."""
    if p == (0,):
        return -1
    return len(p) - 1

def poly_to_str(p):
    """Converts a polynomial to a string representation."""
    if p == (0,): return "0"
    if p == (1,): return "1"
    terms = []
    for i in range(len(p) - 1, -1, -1):
        if p[i] == 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(terms)

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit by searching through polynomials
    of increasing degree.
    """
    x_plus_1 = (1, 1)
    x_4 = (0, 0, 0, 0, 1)
    one = (1,)

    for deg_u in range(1, 10):
        # Iterate through possible degrees for A and B
        max_deg_A = deg_u - 1
        max_deg_B = deg_u

        for dA in range(-1, max_deg_A + 1):
            for dB in range(-1, max_deg_B + 1):
                # We need max(deg(A)+1, deg(B)) == deg_u
                if max(dA + 1, dB) != deg_u:
                    continue

                # Iterate through all polynomials A of degree dA
                coeffs_A_iter = itertools.product([0, 1], repeat=dA + 1) if dA >= 0 else [[0]]
                for coeffs_A_rev in coeffs_A_iter:
                    A = tuple(reversed(coeffs_A_rev))
                    if poly_deg(A) != dA: continue

                    # Iterate through all polynomials B of degree dB
                    coeffs_B_iter = itertools.product([0, 1], repeat=dB + 1) if dB >= 0 else [[0]]
                    for coeffs_B_rev in coeffs_B_iter:
                        B = tuple(reversed(coeffs_B_rev))
                        if poly_deg(B) != dB: continue
                        
                        # Skip trivial unit u=1 (A=0, B=1)
                        if A == (0,) and B == (1,):
                            continue

                        # Check norm equation: A^2(x+1) + A*B*x^4 + B^2 = 1
                        A2 = poly_mul(A, A)
                        B2 = poly_mul(B, B)
                        term1 = poly_mul(A2, x_plus_1)
                        AB = poly_mul(A, B)
                        term2 = poly_mul(AB, x_4)
                        
                        norm = poly_add(poly_add(term1, term2), B2)
                        
                        if norm == one:
                            print(f"The least degree of a non-trivial unit is {deg_u}.")
                            print(f"A unit with this degree is u = A(x)y + B(x), where:")
                            print(f"A(x) = {poly_to_str(A)}")
                            print(f"B(x) = {poly_to_str(B)}")
                            return

find_least_degree_unit()
>>>3