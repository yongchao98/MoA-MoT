import sympy
from sympy import sqrt

def solve_optimization_problem():
    """
    Solves the constrained optimization problem to find the maximum of |b| + |c|
    using the sympy library.
    """
    # We want to maximize b+c subject to b^2 + 4c^2 - 4c = 0 for b,c >= 0.
    # Set up the Lagrange multiplier problem.
    b, c, l = sympy.symbols('b c lambda')

    # Function to maximize
    f = b + c

    # Constraint
    g = b**2 + 4*c**2 - 4*c

    # Gradient equations
    grad_f = [sympy.diff(f, b), sympy.diff(f, c)]
    grad_g = [sympy.diff(g, b), sympy.diff(g, c)]

    # Lagrange equations: grad(f) = lambda * grad(g)
    eq1 = sympy.Eq(grad_f[0], l * grad_g[0])
    eq2 = sympy.Eq(grad_f[1], l * grad_g[1])
    eq3 = sympy.Eq(g, 0)

    # Solve the system of equations
    # We are looking for a solution where b > 0 and c > 0.
    solutions = sympy.solve([eq1, eq2, eq3], (b, c, l), dict=True)
    
    # Filter for the relevant solution (b>0, c>0)
    final_b = 0
    final_c = 0
    max_val = 0
    for sol in solutions:
        # We need to evaluate the numerical values to check conditions
        b_val = sol[b].evalf()
        c_val = sol[c].evalf()
        
        if b_val > 0 and c_val > 0:
            current_max = (sol[b] + sol[c])
            if current_max.evalf() > max_val:
                max_val = current_max
                final_b = sol[b]
                final_c = sol[c]

    # By symmetry, the maximum value for |b|+|c| is the value b+c we found.
    # We can construct one of the extremal polynomials.
    # Let's choose the case b>0, c<0.
    # a = |c|, b = b, c = -|c|
    
    # The optimal values for |b| and |c| are final_b and final_c
    opt_abs_b = final_b
    opt_abs_c = final_c
    
    # We choose one specific case for the signs of a, b, c
    # Let b > 0 and c < 0.
    # The associated 'a' can be a = |c|
    a_val = opt_abs_c
    b_val = opt_abs_b
    c_val = -opt_abs_c
    
    max_sum = opt_abs_b + opt_abs_c

    print("An optimal polynomial is f(x) = a*x^2 + b*x + c, where the coefficients are:")
    print(f"a = {a_val.evalf()} which is {a_val}")
    print(f"b = {b_val.evalf()} which is {b_val}")
    print(f"c = {c_val.evalf()} which is {c_val}")
    print("\nThe maximum value of |b| + |c| is:")
    print(f"{max_sum.evalf()} which is {max_sum}")


solve_optimization_problem()