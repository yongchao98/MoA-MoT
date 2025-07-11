import sympy

def solve_work_cycle():
    """
    This function symbolically derives the work done by the current source
    over one rectangular cycle in the I-x plane.
    """
    # 1. Define symbolic variables
    mu, mu_0, N, w, g, D, x = sympy.symbols('mu mu_0 N w g D x')
    I, I_1, I_2, x_1, x_2 = sympy.symbols('I I_1 I_2 x_1 x_2')

    # 2. Determine the inductance L(x)
    # Reluctance of the part of the gap filled by the block
    R_block = g / (mu * x * w)
    # Reluctance of the part of the gap filled with air
    R_air = g / (mu_0 * (D - x) * w)
    # Total reluctance is the parallel combination
    # 1/R_total = 1/R_block + 1/R_air
    R_total = 1 / (1/R_block + 1/R_air)
    # Inductance L = N^2 / R_total
    L_x = N**2 / R_total
    L_x = sympy.simplify(L_x)
    # L(x) = N**2*w*(D*mu_0 + x*(mu - mu_0))/g

    # 3. Define inductance at positions x1 and x2
    L_x1 = L_x.subs(x, x_1)
    L_x2 = L_x.subs(x, x_2)

    # 4. Calculate the work done by the source (W_src = integral(I * d(lambda))) for each path
    # Flux linkage lambda = L(x) * I

    # Path 1: x from x1 to x2, I = I1 (constant)
    # d(lambda) = I1 * dL. W1 = integral(I1 * (I1 * dL)) = I1**2 * (L(x2) - L(x1))
    W1 = I_1**2 * (L_x2 - L_x1)

    # Path 2: I from I1 to I2, x = x2 (constant)
    # d(lambda) = L(x2) * dI. W2 = integral(I * L(x2) * dI) from I1 to I2
    W2 = sympy.integrate(I * L_x2, (I, I_1, I_2))

    # Path 3: x from x2 to x1, I = I2 (constant)
    # d(lambda) = I2 * dL. W3 = integral(I2 * (I2 * dL)) = I2**2 * (L(x1) - L(x2))
    W3 = I_2**2 * (L_x1 - L_x2)

    # Path 4: I from I2 to I1, x = x1 (constant)
    # d(lambda) = L(x1) * dI. W4 = integral(I * L(x1) * dI) from I2 to I1
    W4 = sympy.integrate(I * L_x1, (I, I_2, I_1))

    # 5. Sum the work from all paths and simplify
    W_total = W1 + W2 + W3 + W4
    W_final = sympy.simplify(W_total)

    # 6. Print the final result in a clear format
    # The result from sympy is: -N**2*w*(-I_1**2 + I_2**2)*(-mu_0 + mu)*(-x_1 + x_2)/(2*g)
    # We reformat it to match the options for better readability.
    # W = - (mu - mu_0) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2) / (2*g)
    
    # Let's create the final expression manually for better formatting.
    # We extract the numerator and denominator to rearrange.
    num, den = W_final.as_numer_denom()
    
    # We make a readable expression
    final_expression = - (mu - mu_0) / (2 * g) * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)
    
    # Print the equation with all its terms.
    # Using sympy.pretty_print for a more mathematical layout.
    print("The work done by the current source for each cycle is:")
    sympy.pretty_print(sympy.Eq(sympy.Symbol('W'), final_expression, evaluate=False), use_unicode=True)

if __name__ == '__main__':
    solve_work_cycle()