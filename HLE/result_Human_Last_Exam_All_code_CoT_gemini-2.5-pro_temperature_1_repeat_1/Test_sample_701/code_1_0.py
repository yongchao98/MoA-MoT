import itertools

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    n1, n2 = len(p1), len(p2)
    n = max(n1, n2)
    p1_padded = p1 + [0] * (n - n1)
    p2_padded = p2 + [0] * (n - n2)
    return [(c1 + c2) % 2 for c1, c2 in zip(p1_padded, p2_padded)]

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if not p1 or not p2:
        return [0]
    n1, n2 = len(p1), len(p2)
    prod = [0] * (n1 + n2 - 1)
    for i in range(n1):
        for j in range(n2):
            prod[i+j] = (prod[i+j] + p1[i] * p2[j]) % 2
    return prod

def poly_deg(p):
    """Calculates the degree of a polynomial."""
    for i in range(len(p) - 1, -1, -1):
        if p[i] == 1:
            return i
    return -1

def poly_to_str(p):
    """Converts a polynomial to a string representation."""
    d = poly_deg(p)
    if d == -1:
        return "0"
    if d == 0:
        return "1"
    terms = []
    for i in range(d, -1, -1):
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
    Finds the unit of least degree by searching systematically.
    """
    x_poly = [0, 1]
    x_plus_1 = [1, 1]
    x_4 = [0, 0, 0, 0, 1]
    one = [1]
    
    for d_unit in range(1, 10):
        print(f"Searching for units of degree {d_unit}...")
        # max(deg(a), deg(b)+1) = d_unit
        max_da = d_unit
        max_db = d_unit - 1

        for da in range(-1, max_da + 1):
            for db in range(-1, max_db + 1):
                if max(da, db + 1) != d_unit:
                    continue

                # Iterate over polynomials a(x) of degree da
                a_coeffs_iter = itertools.product([0, 1], repeat=da + 1) if da >= 0 else [[0]]
                for a_coeffs_tuple in a_coeffs_iter:
                    a = list(a_coeffs_tuple)
                    if poly_deg(a) != da:
                        continue
                    
                    # Iterate over polynomials b(x) of degree db
                    b_coeffs_iter = itertools.product([0, 1], repeat=db + 1) if db >= 0 else [[0]]
                    for b_coeffs_tuple in b_coeffs_iter:
                        b = list(b_coeffs_tuple)
                        if poly_deg(b) != db:
                            continue
                        
                        # Skip trivial unit u=1, where a=1 and b=0
                        if poly_deg(a) == 0 and a[0] == 1 and poly_deg(b) == -1:
                            continue
                        
                        # Norm(u) = a^2 + a*b*x^4 + b^2*(x+1)
                        a_sq = poly_mul(a, a)
                        ab = poly_mul(a, b)
                        abx4 = poly_mul(ab, x_4)
                        b_sq = poly_mul(b, b)
                        b_sq_x1 = poly_mul(b_sq, x_plus_1)
                        
                        norm = poly_add(poly_add(a_sq, abx4), b_sq_x1)
                        
                        if poly_deg(norm) == 0 and norm[0] == 1:
                            print(f"\nFound the unit of least degree: {d_unit}")
                            print(f"The unit is u = a(x) + b(x) * y, where:")
                            print(f"a(x) = {poly_to_str(a)}")
                            print(f"b(x) = {poly_to_str(b)}")
                            return d_unit
    return None

least_degree = find_least_degree_unit()
# The question asks for the least degree, so we would print the number.
# print(f"\nThe least degree of a unit u != 1 is {least_degree}.")