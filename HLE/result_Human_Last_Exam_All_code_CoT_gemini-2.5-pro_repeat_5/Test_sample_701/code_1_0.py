import sympy

def solve():
    """
    Finds the least degree of a unit u != 1 in the ring R by systematically
    searching for units of increasing degree.
    An element u = A(x)*y + B(x) is a unit if its norm is 1, where
    N(u) = A(x)^2*(x+1) + B(x)^2 + A(x)*B(x)*x^4.
    The degree of u is max(deg(A) + 1, deg(B)).
    The script searches for the smallest degree d >= 1 for which a unit exists.
    """
    x = sympy.Symbol('x')
    F2 = sympy.GF(2)
    poly_ring = sympy.Poly(0, x, domain=F2).parent
    x_poly = poly_ring.from_list([1, 0])
    one_poly = poly_ring.from_list([1])
    x_plus_1 = x_poly + one_poly
    x_pow_4 = x_poly**4

    def check_norm(A, B):
        """Computes the norm of u = Ay + B and checks if it is 1."""
        norm = (A * A) * x_plus_1 + B * B + A * B * x_pow_4
        return norm == one_poly

    def get_poly_from_int(n, max_deg):
        """Generates a polynomial in F_2[x] from an integer."""
        if max_deg < 0:
            return poly_ring.zero
        coeffs = []
        temp_n = n
        for _ in range(max_deg + 1):
            coeffs.append(temp_n % 2)
            temp_n //= 2
        coeffs.reverse()
        return poly_ring.from_list(coeffs)

    # Search for the least degree d
    for d in range(1, 6):
        # For a unit of degree d, deg(A) <= d-1 and deg(B) <= d
        max_deg_A = d - 1
        max_deg_B = d
        
        num_A = 2**(max_deg_A + 1) if max_deg_A >= 0 else 1
        num_B = 2**(max_deg_B + 1) if max_deg_B >= 0 else 1

        for i in range(num_A):
            A = get_poly_from_int(i, max_deg_A)
            for j in range(num_B):
                B = get_poly_from_int(j, max_deg_B)
                
                # Skip u=1, which corresponds to A=0, B=1
                if A.is_zero and B == one_poly:
                    continue

                # Check if the degree of u is exactly d
                deg_A = sympy.degree(A) if not A.is_zero else -1
                deg_B = sympy.degree(B) if not B.is_zero else -1
                current_deg = max(deg_A + 1, deg_B)

                if current_deg != d:
                    continue

                # Check if this element is a unit
                if check_norm(A, B):
                    print(d)
                    return

solve()