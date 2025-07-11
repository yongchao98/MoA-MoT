import sympy

def solve_problem():
    """
    This function solves the problem step-by-step using symbolic math.
    """
    # Part 1: Determine the function f(x)
    x, a_coeff, b_coeff, c_coeff = sympy.symbols('x a_coeff b_coeff c_coeff', real=True)
    f_expr = a_coeff * sympy.exp(2*x) + b_coeff * sympy.exp(x) + c_coeff

    # Condition 1: Limit
    # lim_{x->-oo} (f(x)+3)/e^x = 1
    # We solve this manually as explained in the text to get b_coeff=1, c_coeff=-3.
    b_val = 1
    c_val = -3
    f_expr = f_expr.subs([(b_coeff, b_val), (c_coeff, c_val)])
    
    # Condition 2: f(ln(2)) = 0
    # Solve for the coefficient a
    eq1 = f_expr.subs(x, sympy.ln(2))
    a_val = sympy.solve(eq1, a_coeff)[0]
    
    # Final function f(x)
    f = f_expr.subs(a_coeff, a_val)
    print(f"The function is f(x) = {f}")

    # Part 2 & 3: Solve for the integral parameters a and b
    # The integral relation reduces to a = f(ln(b))
    # Let's call the unknowns a_int and b_int to distinguish them.
    a_int, b_int = sympy.symbols('a_int b_int', real=True)
    
    # This is the relation between a and b
    relation = sympy.Eq(a_int, f.subs(x, sympy.ln(b_int)))
    print(f"\nThe integral identity leads to the relation: {relation}")
    
    # To find a unique solution, we assume the implicit condition a=b
    assumption = sympy.Eq(a_int, b_int)
    print(f"Assuming an implicit condition {assumption} to find a unique solution.")

    # Solve the system of equations
    solution = sympy.solve([relation, assumption], (a_int, b_int))

    # The solutions are pairs (a,b). We need b>0 for ln(b) to be real.
    final_a = None
    final_b = None
    for sol_a, sol_b in solution:
        if sol_b > 0:
            final_a = sol_a
            final_b = sol_b
            break

    # Part 4: Final Calculation
    print(f"\nThe solution for the integral parameters is (a, b) = ({final_a}, {final_b})")
    
    a_plus_b = final_a + final_b
    print(f"\nThe value of a + b is: {a_plus_b}")

    # Outputting the numbers in the final equation as requested.
    print("\n--- Final Equation Breakdown ---")
    print(f"f(x) = {sympy.latex(a_val)}*e^(2x) + {b_val}*e^x + {c_val}")
    print(f"The equation to solve is:")
    print(f"integrate(g(x), (x, 0, {final_a})) + integrate(f(x), (x, ln(2), ln({final_b}))) = {final_a}*ln({final_b})")

solve_problem()