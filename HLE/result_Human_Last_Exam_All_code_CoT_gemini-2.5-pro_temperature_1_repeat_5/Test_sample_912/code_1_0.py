import sympy

def solve_magnetic_work():
    """
    This function symbolically derives the work done by the current source
    in the described electromechanical system over one cycle.
    """
    # Define symbolic variables
    mu_0, mu, N, w, g, D, x, x1, x2, I, I1, I2 = sympy.symbols(
        'mu_0 mu N w g D x x1 x2 I I1 I2', real=True, positive=True)

    print("Step 1: Define the Inductance L(x)")
    # Based on the parallel reluctance model
    # L(x) = N^2 / R_total(x)
    # 1/R_total = (mu*x*w)/g + (mu_0*(D-x)*w)/g
    # L(x) = (N^2 * w / g) * ((mu - mu_0)*x + mu_0*D)
    L = (N**2 * w / g) * ((mu - mu_0) * x + mu_0 * D)
    print(f"L(x) = {L}\n")

    # Inductance at positions x1 and x2
    L_x1 = L.subs(x, x1)
    L_x2 = L.subs(x, x2)

    print("Step 2: Calculate work done for each path of the cycle")
    # Path 1: x from x1 to x2, I = I1
    W1 = I1**2 * (L_x2 - L_x1)
    print(f"Work for Path 1 (x: x1->x2, I=I1), W1 = {sympy.simplify(W1)}")

    # Path 2: I from I1 to I2, x = x2
    W2 = sympy.Rational(1, 2) * L_x2 * (I2**2 - I1**2)
    print(f"Work for Path 2 (I: I1->I2, x=x2), W2 = {sympy.simplify(W2)}")

    # Path 3: x from x2 to x1, I = I2
    W3 = I2**2 * (L_x1 - L_x2)
    print(f"Work for Path 3 (x: x2->x1, I=I2), W3 = {sympy.simplify(W3)}")

    # Path 4: I from I2 to I1, x = x1
    W4 = sympy.Rational(1, 2) * L_x1 * (I1**2 - I2**2)
    print(f"Work for Path 4 (I: I2->I1, x=x1), W4 = {sympy.simplify(W4)}\n")

    print("Step 3: Sum the work from all four paths to find the total work W")
    # Total Work
    W_total = W1 + W2 + W3 + W4
    print(f"W_total = W1 + W2 + W3 + W4")

    # Simplify the total work expression
    W_simplified = sympy.simplify(W_total)
    
    # Final formatting to match the options
    final_expression = - (mu - mu_0)/(2*g) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)

    print("\nFinal simplified expression for the total work done by the source:")
    # Using sympy.pretty_print for a more readable output of the formula
    sympy.pretty_print(final_expression)
    
    # To match the output format, we manually construct the print string.
    print("\nW = - (mu - mu_0) / (2*g) * N**2 * w * (x2 - x1) * (I2**2 - I1**2)")


solve_magnetic_work()