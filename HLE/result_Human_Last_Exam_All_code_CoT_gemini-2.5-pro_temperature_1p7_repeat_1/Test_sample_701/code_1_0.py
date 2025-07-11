import itertools

# In Python, we can represent polynomials as lists of coefficients.
# For example, [1, 1, 0, 1] represents 1 + x + x^3 in F_2[x].
# We will use a polynomial class for easier operations.

def poly_to_str(p):
    """Converts a polynomial list to a string."""
    if not p: return "0"
    terms = []
    for i, c in enumerate(p):
        if c == 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(reversed(terms)) if terms else "0"

def poly_add(p1, p2):
    """Adds two polynomials in F_2[x]."""
    res = [0] * max(len(p1), len(p2))
    for i in range(len(p1)):
        res[i] = p1[i]
    for i in range(len(p2)):
        res[i] = (res[i] + p2[i]) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    if not res: return [0]
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials in F_2[x]."""
    if (len(p1)==1 and p1[0]==0) or (len(p2)==1 and p2[0]==0):
        return [0]
    res = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            if p1[i] == 1 and p2[i] == 1:
                res[i+j] = (res[i+j] + p1[i] * p2[j]) % 2
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_deg(p):
    """Returns the degree of a polynomial."""
    if p == [0]: return -1
    return len(p) - 1

def find_unit(max_deg):
    """
    Searches for the unit of least degree up to max_deg.
    A unit u = A(x) + B(x)y has norm 1, where
    N(u) = A^2 + A*B*x^4 + B^2*(x+1) = 1
    """
    x = [0, 1]
    x_plus_1 = [1, 1]
    x_4 = [0, 0, 0, 0, 1]
    one = [1]
    
    for deg in range(1, max_deg + 1):
        num_poly_at_deg = 2**deg
        
        # Iterate over all polynomials A, B with max degree up to `deg`
        # where at least one has degree exactly `deg`.
        for i in range(2**(deg+1)): # A
            A_coeffs = [(i >> k) & 1 for k in range(deg + 1)]
            while len(A_coeffs) > 1 and A_coeffs[-1] == 0: A_coeffs.pop()

            for j in range(2**(deg+1)): # B
                B_coeffs = [(j >> k) & 1 for k in range(deg + 1)]
                while len(B_coeffs) > 1 and B_coeffs[-1] == 0: B_coeffs.pop()
                
                # We want the max degree to be exactly `deg`
                if poly_deg(A_coeffs) != deg and poly_deg(B_coeffs) != deg:
                    continue
                # Skip trivial unit
                if A_coeffs == [1] and B_coeffs == [0]:
                    continue

                A_sq = poly_mul(A_coeffs, A_coeffs)
                B_sq = poly_mul(B_coeffs, B_coeffs)
                AB = poly_mul(A_coeffs, B_coeffs)
                
                term1 = A_sq
                term2 = poly_mul(AB, x_4)
                term3 = poly_mul(B_sq, x_plus_1)
                
                norm = poly_add(poly_add(term1, term2), term3)
                
                if norm == one:
                    print(f"Found a unit of degree {deg}:")
                    print(f"u = ({poly_to_str(A_coeffs)}) + ({poly_to_str(B_coeffs)})y")
                    print("Its norm is 1.")
                    return deg
    return None

min_degree = find_unit(5)
if min_degree:
    print(f"\nThe least degree of a non-trivial unit is {min_degree}.")
else:
    print("No unit found in the specified degree range.")
