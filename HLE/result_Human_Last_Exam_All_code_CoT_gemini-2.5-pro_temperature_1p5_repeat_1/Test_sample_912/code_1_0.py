import sympy

def solve_work_cycle():
    """
    This function symbolically calculates the work done by the current source
    in a magnetic actuator over a rectangular cycle in the I-x plane.
    """
    # 1. Define symbolic variables
    mu_0, mu, N, w, g, D, x = sympy.symbols('mu_0 mu N w g D x', real=True, positive=True)
    I, I_1, I_2 = sympy.symbols('I I_1 I_2', real=True)
    x_1, x_2 = sympy.symbols('x_1 x_2', real=True)

    # 2. Determine the Inductance L(x)
    # The air gap has two parts in parallel.
    # Reluctance of the part with the magnetic block: R_block = g / (mu * x * w)
    # Reluctance of the part with air: R_air = g / (mu_0 * (D - x) * w)
    # Total reluctance of the gap (and the circuit): R_total = 1 / (1/R_block + 1/R_air)
    # The inductance is L(x) = N^2 / R_total.
    L_x = (N**2 * w / g) * (mu_0 * D + (mu - mu_0) * x)
    
    # 3. Define the differential work expression dW = I*L dI + I^2 * (dL/dx) dx
    dL_dx = sympy.diff(L_x, x)
    
    # 4. Integrate over the four paths of the cycle
    
    # Path 1: I = I_1, x from x_1 to x_2. dI = 0.
    # dW1 = I_1^2 * dL/dx * dx
    W1 = sympy.integrate((I_1**2 * dL_dx), (x, x_1, x_2))

    # Path 2: x = x_2, I from I_1 to I_2. dx = 0.
    # dW2 = I * L(x_2) * dI
    L_at_x2 = L_x.subs(x, x_2)
    W2 = sympy.integrate((I * L_at_x2), (I, I_1, I_2))

    # Path 3: I = I_2, x from x_2 to x_1. dI = 0.
    # dW3 = I_2^2 * dL/dx * dx
    W3 = sympy.integrate((I_2**2 * dL_dx), (x, x_2, x_1))

    # Path 4: x = x_1, I from I_2 to I_1. dx = 0.
    # dW4 = I * L(x_1) * dI
    L_at_x1 = L_x.subs(x, x_1)
    W4 = sympy.integrate((I * L_at_x1), (I, I_2, I_1))
    
    # 5. Sum the work from the four paths and simplify
    W_total = sympy.simplify(W1 + W2 + W3 + W4)

    # Print the result in a formatted way.
    # The result from sympy is: -N**2*w*(-I_1**2 + I_2**2)*(-mu_0 + mu)*(-x_1 + x_2)/(2*g)
    # This simplifies to: - (mu - mu_0)/ (2*g) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)
    
    print("The total work done by the current source per cycle is:")
    
    # Explicitly print the formula showing all coefficients mentioned in the final expression
    term1 = "-(mu - mu_0)"
    term2 = "N**2"
    term3 = "w"
    term4 = "(x_2 - x_1)"
    term5 = "(I_2**2 - I_1**2)"
    denominator = "2*g"
    
    print(f"W = ( {term1} * {term2} * {term3} * {term4} * {term5} ) / ( {denominator} )")
    
    # This formula contains the coefficient -1 and the divisor 2.
    print("\nIn standard notation:")
    print(f"W = - ((\u03BC - \u03BC\u2080) * N\u00B2 * w * (x\u2082 - x\u2081) * (I\u2082\u00B2 - I\u2081\u00B2)) / (2g)")


solve_work_cycle()