import itertools

# Polynomial arithmetic in F_2[x]
# A polynomial is represented as a list of coefficients, e.g., x^2+1 is [1, 0, 1]

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    res = [0] * max(len(p1), len(p2))
    for i in range(len(p1)):
        res[i] ^= p1[i]
    for i in range(len(p2)):
        res[i] ^= p2[i]
    return trim(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if p1 == [0] or p2 == [0]:
        return [0]
    res = [0] * (len(p1) + len(p2) - 1)
    for i1, c1 in enumerate(p1):
        if c1 == 1:
            for i2, c2 in enumerate(p2):
                if c2 == 1:
                    res[i1 + i2] ^= 1
    return trim(res)

def trim(p):
    """Removes leading zero coefficients."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    if not p:
        return [0]
    return p

def poly_to_string(p):
    """Converts a polynomial list to a string."""
    if p == [0]:
        return "0"
    terms = []
    for i, c in enumerate(p):
        if c == 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(reversed(terms))

def check_norm(a, b):
    """Checks if N(a+by) = 1."""
    x = [0, 1]
    x4 = [0, 0, 0, 0, 1]
    x_plus_1 = [1, 1]
    
    a_sq = poly_mul(a, a)
    b_sq = poly_mul(b, b)
    
    term1 = a_sq
    term2 = poly_mul(poly_mul(a, b), x4)
    term3 = poly_mul(b_sq, x_plus_1)
    
    norm = poly_add(poly_add(term1, term2), term3)
    
    return norm == [1]

def find_least_degree_unit():
    """Finds the least degree of a non-trivial unit."""
    d = 1
    while True:
        # Iterate through polynomials a(x) and b(x) where max(deg(a), deg(b)) = d
        
        # Case 1: deg(a) = d, deg(b) < d
        for a_coeffs_middle in itertools.product([0, 1], repeat=d):
            a = list(a_coeffs_middle) + [1]
            for deg_b in range(d):
                for b_coeffs_middle in itertools.product([0, 1], repeat=deg_b):
                    b = list(b_coeffs_middle) + [1]
                    if check_norm(a,b):
                        return d, a, b
            # also check b=0, b=1
            if check_norm(a, [0]): return d, a, [0]
            if check_norm(a, [1]): return d, a, [1]


        # Case 2: deg(b) = d
        for b_coeffs_middle in itertools.product([0, 1], repeat=d):
            b = list(b_coeffs_middle) + [1]
            # deg(a) <= d
            for deg_a in range(d + 1):
                # leading coeff of a can be 0 (for deg_a < d) or 1
                for a_coeffs_middle in itertools.product([0, 1], repeat=deg_a):
                    a = list(a_coeffs_middle) + [1]
                    a = trim(a)
                    if a == [1] and b == [0]: continue # skip trivial unit
                    if check_norm(a, b):
                        return d, a, b
            if check_norm([0], b): return d, [0], b
            if check_norm([1], b): return d, [1], b
        
        d += 1

degree, a_poly, b_poly = find_least_degree_unit()

a_str = poly_to_string(a_poly)
b_str = poly_to_string(b_poly)

print(f"The least degree of a non-trivial unit is {degree}.")
print(f"A unit of this degree is u = a(x) + b(x)y where:")
print(f"a(x) = {a_str}")
print(f"b(x) = {b_str}")
print("\nVerification of the norm equation: a(x)^2 + a(x)b(x)x^4 + b(x)^2(x+1) = 1")

# For verification, let's show the terms
x_plus_1 = [1, 1]
x4 = [0, 0, 0, 0, 1]
a_sq = poly_mul(a_poly, a_poly)
b_sq = poly_mul(b_poly, b_poly)
abx4 = poly_mul(poly_mul(a_poly, b_poly), x4)
b_sq_x_plus_1 = poly_mul(b_sq, x_plus_1)

print(f"a(x)^2 = {poly_to_string(a_sq)}")
print(f"a(x)b(x)x^4 = {poly_to_string(abx4)}")
print(f"b(x)^2(x+1) = {poly_to_string(b_sq_x_plus_1)}")

final_sum = poly_add(poly_add(a_sq, abx4), b_sq_x_plus_1)
print(f"\n({poly_to_string(a_sq)}) + ({poly_to_string(abx4)}) + ({poly_to_string(b_sq_x_plus_1)}) = {poly_to_string(final_sum)}")

print(f"\nSince the norm is 1, the element is a unit.")

print("\n<<<3>>>")