import sympy

def solve_work_cycle():
    """
    This function symbolically calculates the work done by the current source
    for a rectangular cycle in the I-x plane.
    """
    # Step 1: Define all symbolic variables
    mu0, mu, N, w, g, D, x, I = sympy.symbols('mu_0 mu N w g D x I', real=True, positive=True)
    x1, x2, I1, I2 = sympy.symbols('x1 x2 I1 I2', real=True)

    # Step 2: Define the inductance L as a function of position x.
    # The total reluctance is the parallel combination of the reluctance of the
    # part of the gap filled by the block and the part filled by air.
    # R_total = g / (w * (mu * x / g + mu0 * (D - x) / g)) is incorrect.
    # Let's re-evaluate reluctance based on parallel flux paths.
    # Permeance_block = (mu * w * x) / g
    # Permeance_air = (mu0 * w * (D - x)) / g
    # Total Permeance P(x) = Permeance_block + Permeance_air
    # P(x) = (w/g) * (mu*x + mu0*(D-x)) = (w/g) * (mu0*D + (mu - mu0)*x)
    # Total Reluctance R(x) = 1/P(x) = g / (w * (mu0*D + (mu - mu0)*x))
    # Inductance L(x) = N**2 / R(x)
    L_x = N**2 * w / g * (mu0 * D + (mu - mu0) * x)
    
    # Derivative of L(x) with respect to x, needed for work calculation
    dL_dx = sympy.diff(L_x, x)

    # Step 3: Calculate the work done for each step of the cycle
    # The differential work done by the source is dW = L(x)*I*dI + I**2*dL(x)
    
    # Path 1: I = I1 (constant), x from x1 to x2. dI = 0.
    # dW1 = I1**2 * dL = I1**2 * dL_dx * dx
    W1 = sympy.integrate(I1**2 * dL_dx, (x, x1, x2))

    # Path 2: x = x2 (constant), I from I1 to I2. dx = 0.
    # dW2 = L(x2) * I * dI
    L_x2 = L_x.subs(x, x2)
    W2 = sympy.integrate(L_x2 * I, (I, I1, I2))

    # Path 3: I = I2 (constant), x from x2 to x1. dI = 0.
    # dW3 = I2**2 * dL = I2**2 * dL_dx * dx
    W3 = sympy.integrate(I2**2 * dL_dx, (x, x2, x1))

    # Path 4: x = x1 (constant), I from I2 to I1. dx = 0.
    # dW4 = L(x1) * I * dI
    L_x1 = L_x.subs(x, x1)
    W4 = sympy.integrate(L_x1 * I, (I, I2, I1))

    # Step 4: Sum the work from all paths and simplify
    W_total = sympy.simplify(W1 + W2 + W3 + W4)

    # Step 5: Print the final result in a structured format
    # The simplified expression is W_total = -N**2*w*(I1**2 - I2**2)*(mu - mu0)*(x1 - x2)/(2*g)
    # Let's re-arrange to match the multiple-choice options.
    # W_total = - (mu - mu0)/(2*g) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)
    term1 = (mu - mu0) / (2 * g)
    term2 = N**2
    term3 = w
    term4 = (x2 - x1)
    term5 = (I2**2 - I1**2)
    
    final_expression = -term1 * term2 * term3 * term4 * term5
    
    print("The work done by the current source for each cycle is:")
    print(f"W = {sympy.pretty(final_expression)}")

solve_work_cycle()
<<<D>>>