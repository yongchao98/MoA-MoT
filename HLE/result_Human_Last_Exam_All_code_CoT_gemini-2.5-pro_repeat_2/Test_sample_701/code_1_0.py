import itertools

# Represents a polynomial in F_2[x] as a list of coefficients (LSB first)
# For example, x^2+1 is represented as [1, 0, 1]

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    res = [0] * max(len(p1), len(p2))
    for i in range(len(p1)):
        res[i] = (res[i] + p1[i]) % 2
    for i in range(len(p2)):
        res[i] = (res[i] + p2[i]) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if p1 == [0] or p2 == [0]:
        return [0]
    res = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            if p1[i] == 1 and p2[j] == 1:
                res[i+j] = (res[i+j] + 1) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def get_poly_deg(p):
    """Gets the degree of a polynomial."""
    if p == [0]:
        return -1
    return len(p) - 1

def poly_to_string(p):
    """Converts a polynomial to a string."""
    if p == [0]:
        return "0"
    terms = []
    for i in range(len(p)):
        if p[i] == 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(reversed(terms))

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit by searching systematically.
    """
    # x+1
    x_plus_1 = [1, 1]
    # x^4
    x_4 = [0, 0, 0, 0, 1]
    
    # Iterate through possible total degrees of the unit
    for degree in range(1, 10):
        # Case 1: deg(a)+1 = degree, deg(b) <= degree
        deg_a = degree - 1
        # Iterate through all polynomials a(x) of degree deg_a
        # A poly of degree d has d+1 coeffs, leading one is 1
        num_coeffs_a = deg_a + 1
        for coeffs_a in itertools.product([0, 1], repeat=num_coeffs_a):
            if deg_a > -1 and coeffs_a[-1] == 0: continue # monic
            a = list(coeffs_a)
            if a == [0]: continue

            # Iterate through all polynomials b(x) of degree <= degree
            for deg_b in range(-1, degree + 1):
                num_coeffs_b = deg_b + 1
                if deg_b == -1:
                    num_coeffs_b = 1 # for the zero poly
                for coeffs_b in itertools.product([0, 1], repeat=num_coeffs_b):
                    if deg_b > -1 and coeffs_b[-1] == 0: continue
                    b = list(coeffs_b)
                    
                    # Exclude trivial unit u=1 (a=0, b=1)
                    if a == [0] and b == [1]:
                        continue
                        
                    # Calculate norm: a^2(x+1) + b^2 + abx^4
                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    ab = poly_mul(a, b)
                    
                    term1 = poly_mul(a_sq, x_plus_1)
                    term3 = poly_mul(ab, x_4)
                    
                    norm = poly_add(poly_add(term1, b_sq), term3)
                    
                    if norm == [1]:
                        print(f"Found a unit u = a(x)y + b(x) with:")
                        print(f"a(x) = {poly_to_string(a)}")
                        print(f"b(x) = {poly_to_string(b)}")
                        unit_deg = max(get_poly_deg(a) + 1, get_poly_deg(b))
                        print(f"The degree of this unit is {unit_deg}.")
                        return unit_deg

        # Case 2: deg(b) = degree, deg(a)+1 < degree
        deg_b = degree
        # Iterate through all polynomials b(x) of degree deg_b
        num_coeffs_b = deg_b + 1
        for coeffs_b in itertools.product([0, 1], repeat=num_coeffs_b):
            if deg_b > -1 and coeffs_b[-1] == 0: continue
            b = list(coeffs_b)

            # Iterate through all polynomials a(x) of degree < degree-1
            for deg_a in range(-1, degree - 1):
                num_coeffs_a = deg_a + 1
                if deg_a == -1:
                    num_coeffs_a = 1
                for coeffs_a in itertools.product([0, 1], repeat=num_coeffs_a):
                    if deg_a > -1 and coeffs_a[-1] == 0: continue
                    a = list(coeffs_a)

                    if a == [0] and b == [1]:
                        continue

                    a_sq = poly_mul(a, a)
                    b_sq = poly_mul(b, b)
                    ab = poly_mul(a, b)
                    
                    term1 = poly_mul(a_sq, x_plus_1)
                    term3 = poly_mul(ab, x_4)
                    
                    norm = poly_add(poly_add(term1, b_sq), term3)
                    
                    if norm == [1]:
                        print(f"Found a unit u = a(x)y + b(x) with:")
                        print(f"a(x) = {poly_to_string(a)}")
                        print(f"b(x) = {poly_to_string(b)}")
                        unit_deg = max(get_poly_deg(a) + 1, get_poly_deg(b))
                        print(f"The degree of this unit is {unit_deg}.")
                        return unit_deg

    return None

if __name__ == '__main__':
    least_degree = find_least_degree_unit()
    if least_degree is None:
        print("No unit found in the searched range.")
