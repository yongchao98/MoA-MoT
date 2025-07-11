import sympy

def solve_overhang_problem():
    """
    This function calculates the maximal stable overhang for three identical cubes
    and determines the integers a, b, c for the expression (a + sqrt(b))/(1 + c).
    """
    print("Derivation of the maximal overhang and the integers a, b, c.")
    print("----------------------------------------------------------")
    # Define the maximal stable offset for a 45-degree rotated cube of side length 1.
    d = sympy.sqrt(2) / 2

    # Define symbols for the centers of the cubes.
    x1, x2, x3 = sympy.symbols('x1 x2 x3')

    # Set up the stability equations.
    # We use the symbol 'd_sym' in equations for clarity in printing.
    d_sym = sympy.Symbol('d')
    eq1 = sympy.Eq(x1 - x2, d_sym)
    eq2 = sympy.Eq((x1 + x2)/2 - x3, d_sym)
    eq3 = sympy.Eq(x1 + x2 + x3, 0)
    
    print("1. The stability equations are set up with d = 1/sqrt(2):")
    print(f"   Equation 1 (C1 on C2):         {eq1}")
    print(f"   Equation 2 ({C1,C2} on C3):      {eq2}")
    print(f"   Equation 3 ({C1,C2,C3} on Table): {eq3}")

    # Solve the system of equations for x1, x2, x3 in terms of d.
    solution = sympy.solve([eq1, eq2, eq3], (x1, x2, x3))
    x1_sol = solution[x1]

    print(f"\n2. Solving the system for x1 (center of the top cube) gives:")
    print(f"   x1 = {sympy.pretty(x1_sol)}")

    # The maximum overhang is the rightmost point of the assembly.
    # If all cubes are rotated 45 degrees, the rightmost edge of cube i is at xi + d.
    # The overhang is max(x1+d, x2+d, x3+d). Since x1 is the largest position,
    # the max overhang is determined by the top cube.
    max_overhang_expr = x1_sol + d_sym
    
    print(f"\n3. The maximal overhang is the rightmost edge of the top cube, at x1 + d:")
    print(f"   Maximal Overhang = x1 + d = {sympy.pretty(max_overhang_expr)}")
    
    # Substitute the value of d = 1/sqrt(2) to get the final answer.
    max_overhang_val = max_overhang_expr.subs(d_sym, d)
    
    print(f"\n4. Substituting d = 1/sqrt(2), the overhang is:")
    print(f"   Maximal Overhang = {sympy.pretty(max_overhang_val)}")

    # Express the result in the required format (a + sqrt(b)) / (1 + c).
    # The overhang is (11*sqrt(2))/12 = (0 + sqrt(121*2))/12 = (0 + sqrt(242))/12.
    num, den = max_overhang_val.as_numer_denom()
    
    a = 0
    # num is 11*sqrt(2) = sqrt(121*2) = sqrt(242)
    b = num**2 / 2 * 2 # To handle sympy expression. Simplified: 11**2 * 2
    b = int(b)
    c = int(den) - 1

    final_expr_str = f"({a} + sqrt({b})) / (1 + {c})"
    print(f"\n5. We format this result as (a + sqrt(b)) / (1 + c):")
    print(f"   {sympy.pretty(max_overhang_val)} = {final_expr_str}")
    
    print("\n6. The non-negative integers a, b, c that satisfy the conditions are:")
    print(f"   a = {a}")
    print(f"   b = {b}")
    print(f"   c = {c}")

if __name__ == '__main__':
    solve_overhang_problem()
