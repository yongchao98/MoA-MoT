import sympy

def solve_fourier_restriction_problem():
    """
    This script solves for the largest p such that no L^p function on R^3
    can have its Fourier support on the moment curve.

    The solution relies on the sharp boundedness theorem for the extension
    operator associated with the moment curve in R^3.
    A non-zero function f in L^p with Fourier support on the curve exists
    if and only if the extension operator E: L^q -> L^p is bounded for some q.

    The conditions for boundedness are:
    1) p >= q
    2) 3/p + 2/q <= 2
    """

    p, q = sympy.symbols('p q', positive=True)

    # We want to find the threshold for p. A non-zero function exists if we can
    # find *any* q >= 1 that satisfies the boundedness conditions.
    # We combine the inequalities to find the absolute lower bound on p for
    # boundedness to be possible.

    # From condition (1), p >= q, we have 1/p <= 1/q.
    # Substitute this into condition (2):
    # 3/p + 2/q <= 2
    # Since 1/q >= 1/p, the left side is smallest when q is largest.
    # To find the limit on p, we use the tightest constraint from p >= q,
    # which is when q approaches p.
    # 3/p + 2/p <= 3/p + 2/q <= 2
    # This gives the necessary condition on p:
    # 5/p <= 2  =>  p >= 5/2

    # This means that for boundedness to be possible at all, we must have p >= 2.5.
    # For any p > 2.5, we can choose q = p. This satisfies p >= q.
    # The second condition becomes 3/p + 2/p = 5/p < 2, which is true.
    # So for any p > 2.5, a non-zero function exists.

    # The threshold is therefore at p = 5/2. This is the largest value of p
    # for which no non-zero function can be constructed this way.

    # The final equation for the threshold is:
    final_equation_str = "5 / p = 2"
    print("The problem reduces to finding the threshold for the boundedness of the Fourier extension operator.")
    print("The necessary condition for boundedness derived from the sharp theorem is 5/p <= 2.")
    print("The threshold value of p is found by solving the equation at the boundary.")
    print(f"Final equation: {final_equation_str}")

    # Solve the equation
    p_val = sympy.solve(sympy.Eq(5/p, 2), p)
    
    # Output the numbers in the final equation and the solution
    numerator = 5
    denominator_var = 'p'
    rhs = 2
    
    print(f"Solving for {denominator_var} in the equation: {numerator} / {denominator_var} = {rhs}")
    
    result = p_val[0]
    print(f"{denominator_var} = {numerator} / {rhs}")
    print(f"{denominator_var} = {float(result)}")

solve_fourier_restriction_problem()