import sympy

def solve_work_cycle():
    """
    This function symbolically calculates the work done by the current source
    in a rectangular I-x cycle for the given magnetic actuator.
    """
    # Define the symbolic variables
    mu, mu0, N, w, g, D, x = sympy.symbols('mu mu_0 N w g D x', real=True, positive=True)
    I, x1, x2, I1, I2 = sympy.symbols('I x_1 x_2 I_1 I_2', real=True)

    # Step 1 & 2: Calculate Reluctance and Inductance
    # The air gap has two parallel magnetic paths:
    # Path 1: Through the block, with area x*w, length g, permeability mu
    # Path 2: Through air, with area (D-x)*w, length g, permeability mu0
    # Reluctance R = length / (permeability * area)
    # The total permeance (1/R_total) is the sum of the individual permeances.
    permeance_total = (mu * x * w) / g + (mu0 * (D - x) * w) / g
    
    # Inductance L(x) = N^2 * Permeance
    L_x = N**2 * permeance_total
    
    # Let's simplify the inductance expression
    L_x = sympy.simplify(L_x)
    # L_x = N**2*w*(D*mu0 + x*(mu - mu0))/g

    # Step 3: Define the work integral for each path of the cycle
    # dW = I * d(lambda) = I * d(L*I) = I * (L*dI + I*dL) = L*I*dI + I^2*dL
    # dL = (dL/dx) * dx
    dL_dx = sympy.diff(L_x, x)

    # Step 4: Integrate over the four steps of the cycle
    # Step 1: x from x1 to x2, I = I1 (dI=0)
    W1 = sympy.integrate(I1**2 * dL_dx, (x, x1, x2))

    # Step 2: I from I1 to I2, x = x2 (dx=0)
    L_at_x2 = L_x.subs(x, x2)
    W2 = sympy.integrate(L_at_x2 * I, (I, I1, I2))

    # Step 3: x from x2 to x1, I = I2 (dI=0)
    W3 = sympy.integrate(I2**2 * dL_dx, (x, x2, x1))

    # Step 4: I from I2 to I1, x = x1 (dx=0)
    L_at_x1 = L_x.subs(x, x1)
    W4 = sympy.integrate(L_at_x1 * I, (I, I2, I1))

    # Step 5: Sum the work from all steps
    W_total = W1 + W2 + W3 + W4

    # Simplify the final expression
    W_total_simplified = sympy.simplify(W_total)
    
    # The problem asks to output the numbers in the equation.
    # The symbolic result is W = - (mu - mu_0)/(2*g) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)
    # The numbers are -1, 2, 2, 2, 1, 1. Let's print the formula in a clear way.
    
    print("The derived total work done by the source per cycle is:")
    # We use sympy.pretty_print to display the equation in a more readable format.
    W = sympy.Symbol('W')
    final_equation = sympy.Eq(W, W_total_simplified)
    sympy.pretty_print(final_equation)

solve_work_cycle()