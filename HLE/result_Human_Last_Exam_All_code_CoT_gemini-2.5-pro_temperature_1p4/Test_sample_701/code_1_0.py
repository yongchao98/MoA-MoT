import itertools

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x] represented as lists of coefficients."""
    lmax = max(len(p1), len(p2))
    res = [0] * lmax
    for i in range(lmax):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res[i] = (c1 + c2) % 2
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
            res[i+j] = (res[i+j] + p1[i] * p2[j]) % 2
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_to_str(p):
    """Converts a polynomial list to a string."""
    if p == [0]: return "0"
    if p == [1]: return "1"
    terms = []
    for i in range(len(p)):
        if p[i] == 1:
            if i == 0: terms.append("1")
            elif i == 1: terms.append("x")
            else: terms.append(f"x^{i}")
    return " + ".join(reversed(terms))

def check_unit_norm(a, b):
    """Checks if a(x)^2 + a(x)b(x)x^4 + b(x)^2(x+1) = 1."""
    x_plus_1 = [1, 1]
    x_4 = [0, 0, 0, 0, 1]
    
    a_sq = poly_mul(a, a)
    b_sq = poly_mul(b, b)
    
    term1 = a_sq
    term2 = poly_mul(poly_mul(a, b), x_4)
    term3 = poly_mul(b_sq, x_plus_1)
    
    norm = poly_add(poly_add(term1, term2), term3)
    return norm == [1]

def find_unit():
    """Searches for the unit of least degree."""
    max_total_degree = 1
    while True:
        # degree(u) = max(deg(a), deg(b)+1)
        for d in range(1, max_total_degree + 1):
            # Iterate through all possible degrees for a and b
            for deg_b in range(d):
                deg_a = d if deg_b + 1 < d else d # deg(a) can be up to d
                
                # Check pairs (a,b) with deg(a)=deg_a, deg(b)=deg_b
                # (and other way around if d=deg_b+1)
                
                # Let's just do a simpler, exhaustive search by degree
                pass

        # Simplified Search Strategy
        deg_u = max_total_degree
        # deg_b goes from 0 up to deg_u - 1
        for deg_b in range(deg_u):
            # deg_a goes from 0 up to deg_u
            for deg_a in range(deg_u + 1):
                if max(deg_a, deg_b + 1) != deg_u:
                    continue

                # Generate all polynomials for deg_a and deg_b
                for coeffs_a in itertools.product([0, 1], repeat=deg_a):
                    a = list(coeffs_a) + [1]
                    for coeffs_b in itertools.product([0, 1], repeat=deg_b):
                        b = list(coeffs_b) + [1]
                        
                        if check_unit_norm(a, b):
                           print(f"Found a unit u = a(x) + b(x)y with:")
                           print(f"a(x) = {poly_to_str(a)}")
                           print(f"b(x) = {poly_to_str(b)}")
                           print(f"The degree of this unit is max(deg(a), deg(b)+1) = max({deg_a}, {deg_b}+1) = {deg_u}.")
                           return deg_u
        
        max_total_degree += 1
        if max_total_degree > 5: # Safety break
            return -1


result = find_unit()
print("The least degree of a non-trivial unit is:")
print(result)
