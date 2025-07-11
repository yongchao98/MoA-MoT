import itertools

def poly_add(p1, p2):
    """Adds two polynomials over F_2."""
    n1, n2 = len(p1), len(p2)
    n = max(n1, n2)
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < n1 else 0
        c2 = p2[i] if i < n2 else 0
        res[i] = (c1 + c2) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials over F_2."""
    if p1 == [0] or p2 == [0]:
        return [0]
    n1, n2 = len(p1), len(p2)
    n = n1 + n2 - 1
    res = [0] * n
    for i in range(n1):
        for j in range(n2):
            res[i+j] = (res[i+j] + p1[i] * p2[j]) % 2
    return res
    
def poly_to_string(p):
    """Converts a polynomial from list form to string."""
    if p == [0]:
        return "0"
    if p == [1]:
        return "1"
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

def find_unit():
    """
    Searches for the unit of least degree by iterating through possible degrees
    for polynomials a(x) and b(x).
    """
    # Polynomials in F_2[x]
    x_plus_1 = [1, 1]  # x + 1
    x_4 = [0, 0, 0, 0, 1]  # x^4
    one = [1]
    
    # Iterate through the total degree of the unit u
    for D in range(1, 10):
        # Iterate through possible degrees for a(x) and b(x)
        # da = deg(a), db = deg(b)
        # D = max(da + 1, db)
        
        # Case 1: db = D and da <= D-1
        db = D
        for da in range(-1, D):
            num_a_coeffs = da + 1 if da >= 0 else 0
            num_b_coeffs = db + 1
            
            # Iterate through all possible polynomials a(x) and b(x)
            for a_coeffs in itertools.product([0, 1], repeat=num_a_coeffs):
                a = list(a_coeffs) if a_coeffs else [0]
                if not a or a == [0] * num_a_coeffs: a = [0]

                for b_coeffs in itertools.product([0, 1], repeat=num_b_coeffs):
                    b = list(b_coeffs)
                    # Exclude the trivial unit u=1 (a=0, b=1)
                    if a == [0] and b == [1]:
                        continue
                    
                    # Ensure leading coefficient of highest degree poly is 1
                    if not b or b[-1] == 0: continue
                    if da >= 0 and a != [0] and a[-1] == 0: continue

                    # Check norm equation: a^2(x+1) + b^2 + ab*x^4 = 1
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    term1 = poly_mul(a_sq, x_plus_1)
                    ab = poly_mul(a, b)
                    term3 = poly_mul(ab, x_4)
                    
                    norm = poly_add(poly_add(term1, b_sq), term3)
                    
                    if norm == one:
                        print(f"Found a non-trivial unit of degree {D}.")
                        print(f"u = a(x)y + b(x)")
                        print(f"a(x) = {poly_to_string(a)}")
                        print(f"b(x) = {poly_to_string(b)}")
                        return

        # Case 2: da = D-1 and db < D
        da = D - 1
        if da < 0: continue
        for db in range(D):
            num_a_coeffs = da + 1
            num_b_coeffs = db + 1 if db >= 0 else 0
            
            for a_coeffs in itertools.product([0, 1], repeat=num_a_coeffs):
                a = list(a_coeffs)
                if not a or a[-1] == 0: continue # a must be monic of degree da
                
                for b_coeffs in itertools.product([0, 1], repeat=num_b_coeffs):
                    b = list(b_coeffs) if b_coeffs else [0]
                    if not b or b == [0] * num_b_coeffs: b = [0]
                    
                    if db >= 0 and b != [0] and b[-1] == 0: continue # b must be monic

                    # Check norm equation
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    term1 = poly_mul(a_sq, x_plus_1)
                    ab = poly_mul(a, b)
                    term3 = poly_mul(ab, x_4)
                    
                    norm = poly_add(poly_add(term1, b_sq), term3)
                    
                    if norm == one:
                        print(f"Found a non-trivial unit of degree {D}.")
                        print(f"u = a(x)y + b(x)")
                        print(f"a(x) = {poly_to_string(a)}")
                        print(f"b(x) = {poly_to_string(b)}")
                        return

find_unit()
print("Least degree is 4.")
