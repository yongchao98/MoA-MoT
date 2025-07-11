import sympy

def find_unit_degree():
    """
    Searches for the least degree of a non-trivial unit in the ring
    R = F_2[x,y]/(y^2 + x^4*y + x + 1).

    The norm of an element u = a(x) + b(x)y is N(u) = a^2 + a*b*x^4 + b^2*(x+1).
    A non-trivial unit u != 1 exists if N(u) = 1 for (a,b) != (1,0).
    The degree of u is max(deg(a), deg(b)+1).
    """
    x = sympy.Symbol('x')
    F2 = sympy.GF(2)

    # Convert integer to a polynomial in F_2[x]
    def int_to_poly(n, deg):
        if n == 0:
            return sympy.Poly(0, x, domain=F2)
        coeffs = []
        temp_n = n
        for _ in range(deg + 1):
            coeffs.append(temp_n % 2)
            temp_n //= 2
        return sympy.Poly(coeffs[::-1], x, domain=F2)

    # The constant polynomial 1 in F_2[x]
    one_poly = sympy.Poly(1, x, domain=F2)
    x_poly = sympy.Poly(x, domain=F2)
    x4_poly = sympy.Poly(x**4, domain=F2)

    # Loop through possible degrees for the unit u
    for unit_deg in range(1, 5): # Search up to degree 4
        # Case 1: deg(b) + 1 = unit_deg, deg(a) <= unit_deg
        deg_b = unit_deg - 1
        if deg_b >= 0:
            for b_int in range(1, 2**(deg_b + 1)): # b != 0
                b = int_to_poly(b_int, deg_b)
                if b.degree() != deg_b:
                    continue
                
                max_deg_a = unit_deg
                for a_int in range(2**(max_deg_a + 1)):
                    a = int_to_poly(a_int, max_deg_a)
                    
                    # Norm calculation
                    norm = a**2 + a * b * x4_poly + b**2 * (x_poly + one_poly)
                    
                    if norm == one_poly:
                        print(f"Found a unit of degree {unit_deg}:")
                        print(f"u = ({a.as_expr()}) + ({b.as_expr()})*y")
                        return unit_deg

        # Case 2: deg(a) = unit_deg, deg(b) + 1 < unit_deg
        deg_a = unit_deg
        max_deg_b = unit_deg - 2
        for a_int in range(2**deg_a, 2**(deg_a + 1)):
            a = int_to_poly(a_int, deg_a)
            if a.degree() != deg_a:
                continue

            for b_int in range(2**(max_deg_b + 1)) if max_deg_b >= 0 else [0]:
                 b = int_to_poly(b_int, max_deg_b)

                 # Norm calculation
                 norm = a**2 + a * b * x4_poly + b**2 * (x_poly + one_poly)
                 
                 if norm == one_poly:
                     print(f"Found a unit of degree {unit_deg}:")
                     print(f"u = ({a.as_expr()}) + ({b.as_expr()})*y")
                     return unit_deg
    return None

result_degree = find_unit_degree()
if result_degree:
    print(f"\nThe least degree of a unit u != 1 is {result_degree}.")
