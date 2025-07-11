import sympy as sp

def solve_work_cycle():
    """
    This function symbolically derives the work done by the current source
    over a rectangular cycle in the I-x plane for the given magnetic actuator.
    """
    # Define the symbolic variables for the physical quantities
    mu, mu_0, N, g, w, D, x, I = sp.symbols('mu mu_0 N g w D x I', real=True, positive=True)
    x1, x2, I1, I2 = sp.symbols('x_1 x_2 I_1 I_2', real=True, positive=True)

    # Step 1: Define the magnetic reluctance of the two parallel paths in the gap
    # Reluctance of the path through the magnetic block
    R_block = g / (mu * x * w)
    # Reluctance of the path through the air
    R_air = g / (mu_0 * (D - x) * w)

    # The total reluctance is the parallel combination of the two
    R_total = sp.simplify(1 / (1/R_block + 1/R_air))

    # Step 2: Calculate the inductance L(x)
    L_x = N**2 / R_total
    L_x = sp.simplify(L_x)

    # Step 3: Formulate the work integral components
    # The work done by the source is W = integral(I * d(lambda)), where lambda = L(x)*I
    # d(lambda) = L(x)*dI + I*(dL/dx)*dx
    # W = integral(I*L(x) dI) + integral(I**2 * (dL/dx) dx)

    # Step 4: Evaluate the work for each path of the cycle
    # Path 1: x from x1 to x2, I = I1 (constant, so dI=0)
    dL_dx = sp.diff(L_x, x)
    integrand1 = (I**2 * dL_dx).subs(I, I1)
    W1 = sp.integrate(integrand1, (x, x1, x2))

    # Path 2: I from I1 to I2, x = x2 (constant, so dx=0)
    integrand2 = I * L_x.subs(x, x2)
    W2 = sp.integrate(integrand2, (I, I1, I2))

    # Path 3: x from x2 to x1, I = I2 (constant, so dI=0)
    integrand3 = (I**2 * dL_dx).subs(I, I2)
    W3 = sp.integrate(integrand3, (x, x2, x1))

    # Path 4: I from I2 to I1, x = x1 (constant, so dx=0)
    integrand4 = I * L_x.subs(x, x1)
    W4 = sp.integrate(integrand4, (I, I2, I1))

    # Total work is the sum of the work done on all paths
    W_total = sp.simplify(W1 + W2 + W3 + W4)

    # Step 5: Present the final result
    print("The derived formula for the work done by the current source per cycle is:")
    
    # Deconstruct the expression for clear printing
    # The expression is of the form: - (mu - mu_0) * N**2 * w * (x2 - x1) * (I2**2 - I1**2) / (2*g)
    term_mu = sp.pretty(mu - mu_0)
    term_N_sq = sp.pretty(N**2)
    term_w = sp.pretty(w)
    term_x = sp.pretty(x2 - x1)
    term_I_sq = sp.pretty(I2**2 - I1**2)
    term_g = sp.pretty(g)
    
    print("\nW = - ( ({} * {} * {} * {} * {}) / (2 * {}) )".format(
        term_mu, term_N_sq, term_w, term_x, term_I_sq, term_g
    ))

    # Pretty print the final symbolic formula
    print("\nSymbolic Result:")
    sp.pprint(W_total, use_unicode=True)
    
    # Match with the given options
    # A. W = (mu_0*N**2*w*(x2 - x1)*(I2**2 + I1**2))/(2*g)
    # B. W = - (mu + mu_0)*N**2*w*(x2 - x1)*(I2**2 - I1**2)/(2*g)
    # C. W = (mu - mu_0)*N**2*w*(x2 + x1)*(I2 - I1)/g
    # D. W = - (mu - mu_0)*N**2*w*(x2 - x1)*(I2**2 - I1**2)/(2*g)
    # E. W = - (mu - mu_0)*N**2*w*(x2 + x1)*(I2**2 - I1**2)/(2*g)
    # Our result W_total = -N**2*w*(-I1**2 + I2**2)*(-x1 + x2)*(-mu_0 + mu)/(2*g)
    # This is identical to option D.

if __name__ == '__main__':
    solve_work_cycle()
    print("\nBased on the derivation, the correct option is D.")
    print("<<<D>>>")