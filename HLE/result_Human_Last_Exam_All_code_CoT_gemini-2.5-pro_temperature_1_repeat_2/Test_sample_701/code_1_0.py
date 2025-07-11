import numpy as np

# Use numpy's polynomial class for convenience, with coefficients in F_2
# We will handle coefficients by taking them modulo 2.
def f2_poly(coeffs):
    """Represents a polynomial in F_2[x]."""
    p = np.polynomial.Polynomial(coeffs)
    # Ensure coefficients are always 0 or 1
    p.coef = np.mod(p.coef, 2)
    return p

def poly_to_str(p):
    """Custom string representation for F_2 polynomials."""
    if p.degree() < 0:
        return "0"
    s = []
    for i in range(p.degree(), -1, -1):
        if p.coef[i] == 1:
            if i == 0:
                s.append("1")
            elif i == 1:
                s.append("x")
            else:
                s.append(f"x^{i}")
    return " + ".join(s) if s else "0"

def generate_polynomials(degree):
    """Generate all polynomials of a given degree in F_2[x]."""
    if degree < 0:
        yield f2_poly([])
        return
    if degree == 0:
        yield f2_poly([0])
        yield f2_poly([1])
        return
    
    # Iterate through all possible coefficient combinations
    for i in range(2**degree):
        coeffs = []
        temp_i = i
        for _ in range(degree):
            coeffs.append(temp_i % 2)
            temp_i //= 2
        # Add the leading coefficient, which must be 1
        coeffs.append(1)
        yield f2_poly(coeffs)

def find_least_degree_unit():
    """
    Searches for the least degree unit u != 1 in the ring R.
    """
    x = f2_poly([0, 1])
    x_plus_1 = f2_poly([1, 1])
    x4 = x**4
    
    # Loop through possible degrees for the unit u
    for degree in range(1, 10):
        print(f"Checking for units of degree {degree}...")
        # da = deg(a), db = deg(b)
        # degree = max(da, db + 1)
        
        # Iterate through possible degrees of a and b
        for da in range(degree + 1):
            for db in range(degree):
                if max(da, db + 1) != degree:
                    continue

                # Generate all polynomials for the given degrees
                a_candidates = list(generate_polynomials(da)) if da >= 0 else [f2_poly([0])]
                b_candidates = list(generate_polynomials(db)) if db >= 0 else [f2_poly([0])]

                for a in a_candidates:
                    for b in b_candidates:
                        # A unit must satisfy a(1)=1 and b(1)=0
                        if np.mod(a(1), 2) != 1 or np.mod(b(1), 2) != 0:
                            continue

                        # We are looking for u != 1, so b cannot be the zero polynomial
                        if b.degree() < 0:
                            continue

                        # Calculate the norm: N(u) = a^2 + a*b*x^4 + b^2*(x+1)
                        a_sq = a * a
                        b_sq = b * b
                        term1 = a_sq
                        term2 = a * b * x4
                        term3 = b_sq * x_plus_1
                        
                        norm = term1 + term2 + term3
                        norm.coef = np.mod(norm.coef, 2)
                        
                        # Check if the norm is 1
                        if norm.degree() == 0 and norm.coef[0] == 1:
                            print("\nFound a unit!")
                            print(f"u = a(x) + b(x)y")
                            print(f"a(x) = {poly_to_str(a)}")
                            print(f"b(x) = {poly_to_str(b)}")
                            print(f"Degree of the unit is {degree}.")
                            return degree

    return None

if __name__ == '__main__':
    least_degree = find_least_degree_unit()
    if least_degree is None:
        print("\nCould not find a unit within the searched degrees.")
    else:
        print(f"\n<<< {least_degree} >>>")
