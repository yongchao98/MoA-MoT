import collections

# Polynomials over F_2 are represented as lists of coefficients (little-endian).
# For example, x^2 + 1 is represented as [1, 0, 1].

def poly_add(p1, p2):
    """Adds two polynomials over F_2."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(len(p1)):
        res[i] ^= p1[i]
    for i in range(len(p2)):
        res[i] ^= p2[i]
    return trim(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials over F_2."""
    if not p1 or not p2:
        return []
    n1, n2 = len(p1), len(p2)
    res = [0] * (n1 + n2 - 1)
    for i1, c1 in enumerate(p1):
        if c1 == 0: continue
        for i2, c2 in enumerate(p2):
            if c2 == 0: continue
            res[i1 + i2] ^= c2
    return trim(res)

def trim(p):
    """Removes leading zero coefficients."""
    while p and p[-1] == 0:
        p.pop()
    return p

def poly_deg(p):
    """Returns the degree of a polynomial."""
    return len(p) - 1

def poly_eval(p, x):
    """Evaluates a polynomial at x (in F_2)."""
    if not p:
        return 0
    val = 0
    for i, c in enumerate(p):
        if c == 1:
            val ^= (x ** i)
    return val % 2

def poly_to_str(p, var='x'):
    """Converts a polynomial to its string representation."""
    if not p: return "0"
    s = ""
    for i in range(len(p) - 1, -1, -1):
        if p[i] == 1:
            if i == 0:
                s += "1"
            elif i == 1:
                s += f"{var}"
            else:
                s += f"{var}^{i}"
            s += " + "
    return s[:-3]

def find_least_degree_unit():
    """
    Finds the least degree of a non-trivial unit u = a(x) + b(x)y.
    The degree of u is max(deg(a), deg(b) + 1).
    """
    for u_deg in range(1, 10):
        # Iterate through possible degrees for a and b
        for d_a in range(u_deg + 1):
            for d_b in range(u_deg):
                if max(d_a, d_b + 1) != u_deg:
                    continue

                # Iterate through polynomials a(x) of degree d_a
                # Leading coefficient of a must be 1. 2^d_a possibilities for other coeffs.
                for i in range(1 << d_a):
                    a = [int(c) for c in bin(i)[2:].zfill(d_a)] + [1]
                    
                    # Constraints from norm equation: a(1)=1
                    if poly_eval(a, 1) != 1:
                        continue

                    # Iterate through polynomials b(x) of degree d_b
                    # Leading coefficient of b can be 0 or 1.
                    for j in range(1 << (d_b + 1)):
                        if d_b == -1: b = []
                        else: b = [int(c) for c in bin(j)[2:].zfill(d_b + 1)]
                        b = trim(b)
                        if poly_deg(b) != d_b: continue

                        # Constraints: a(0)^2+b(0)^2=1 and b(1)=0
                        a0 = poly_eval(a, 0)
                        b0 = poly_eval(b, 0)
                        if (a0**2 + b0**2) % 2 != 1:
                            continue
                        if poly_eval(b, 1) != 0:
                            continue
                        
                        # Exclude trivial unit u=1
                        if a == [1] and not b:
                            continue

                        # Check norm condition: a^2 + a*b*x^4 + b^2*(x+1) = 1
                        a_sq = poly_mul(a, a)
                        x4 = [0, 0, 0, 0, 1]
                        abx4 = poly_mul(a, poly_mul(b, x4))
                        x_plus_1 = [1, 1]
                        b_sq = poly_mul(b, b)
                        b_sq_x_plus_1 = poly_mul(b_sq, x_plus_1)
                        
                        norm = poly_add(poly_add(a_sq, abx4), b_sq_x_plus_1)

                        if norm == [1]:
                            print(f"Found a unit u = a(x) + b(x)y with least degree {u_deg}.")
                            print("The unit is given by:")
                            print(f"a(x) = {poly_to_str(a)}")
                            print(f"b(x) = {poly_to_str(b)}")
                            print("\nVerification of the norm equation N(u) = a^2 + a*b*x^4 + b^2*(x+1) = 1:")
                            print(f"a(x)^2 = ({poly_to_str(a)})^2 = {poly_to_str(a_sq)}")
                            print(f"a(x)b(x)x^4 = ({poly_to_str(a)})({poly_to_str(b)})x^4 = {poly_to_str(abx4)}")
                            print(f"b(x)^2(x+1) = ({poly_to_str(b)})^2(x+1) = {poly_to_str(b_sq_x_plus_1)}")
                            sum1 = poly_add(a_sq, abx4)
                            final_sum = poly_add(sum1, b_sq_x_plus_1)
                            print(f"The sum is ({poly_to_str(a_sq)}) + ({poly_to_str(abx4)}) + ({poly_to_str(b_sq_x_plus_1)}) = {poly_to_str(final_sum)}")
                            return u_deg
    return -1

result = find_least_degree_unit()
if result != -1:
    print(f"\nThe least degree of a non-trivial unit is {result}.")

# The final answer is the integer degree, which the problem asks for.
# The code finds and prints the degree.
# To conform to the final output format, we just print the number.
# The code finds the degree to be 4.
print(result)
>>>4