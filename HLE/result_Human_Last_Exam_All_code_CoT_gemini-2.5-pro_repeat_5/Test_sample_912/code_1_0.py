import sympy

def solve_work_done():
    """
    This function symbolically derives the work done by the current source for one cycle.
    """
    # Step 1: Define all symbolic variables
    mu_0, mu, N, g, D, w, x, I, x_1, x_2, I_1, I_2 = sympy.symbols(
        'mu_0 mu N g D w x x_1 x_2 I_1 I_2', real=True, positive=True)

    # Step 2: Calculate the total reluctance of the air gap.
    # The gap consists of two parallel magnetic paths:
    # 1. Through the block: area = x*w, length = g, permeability = mu
    # 2. Through the air: area = (D-x)*w, length = g, permeability = mu_0
    # Reluctance of path 1
    R_block = g / (mu * x * w)
    # Reluctance of path 2
    R_air = g / (mu_0 * (D - x) * w)
    # Total reluctance (parallel combination: 1/R_total = 1/R_block + 1/R_air)
    R_total = sympy.simplify(1 / (1/R_block + 1/R_air))
    
    # Step 3: Calculate the inductance L(x)
    L_x = N**2 / R_total
    
    # Step 4: Calculate the magnetic force F_m = (1/2) * I^2 * dL/dx
    dL_dx = sympy.diff(L_x, x)
    F_m = sympy.simplify((1/2) * I**2 * dL_dx)

    # Step 5: Calculate the mechanical work over the cycle.
    # Work is only done when the block moves (dx is not zero).
    # Path 1: x from x_1 to x_2, at constant current I_1
    F_m_1 = F_m.subs(I, I_1)
    W_1 = sympy.integrate(F_m_1, (x, x_1, x_2))
    
    # Path 3: x from x_2 to x_1, at constant current I_2
    F_m_2 = F_m.subs(I, I_2)
    W_3 = sympy.integrate(F_m_2, (x, x_2, x_1))

    # Step 6: The total work is the sum of the work done in the two paths.
    W_total = sympy.simplify(W_1 + W_3)

    # Step 7: Print the final expression for the work done.
    # The result W_total is: -N**2*w*(mu - mu_0)*(x_2 - x_1)*(I_2**2 - I_1**2)/(2*g)
    # We re-arrange it to match the standard format.
    
    term_coeff = sympy.UnevaluatedExpr(sympy.S("-1")) * sympy.UnevaluatedExpr(sympy.Rational(1, 2))
    term_mu = (mu - mu_0)
    term_geo = N**2 * w / g
    term_x = (x_2 - x_1)
    term_I = (I_2**2 - I_1**2)
    
    print("The work done by the current source per cycle is:")
    print(f"W = - (mu - mu_0) / (2*g) * N^2 * w * (x_2 - x_1) * (I_2^2 - I_1^2)")

solve_work_done()