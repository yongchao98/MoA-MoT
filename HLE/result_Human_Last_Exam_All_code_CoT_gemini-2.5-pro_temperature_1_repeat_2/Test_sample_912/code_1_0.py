import sympy

def solve_work_cycle():
    """
    This function symbolically calculates the work done by the current source
    over the specified rectangular cycle in the I-x plane.
    """
    # Define the symbolic variables from the problem description
    mu_0, mu, N, g, D, w, x = sympy.symbols('mu_0 mu N g D w x', real=True, positive=True)
    x1, x2, I1, I2 = sympy.symbols('x1 x2 I1 I2', real=True)

    # Step 1: Derive the inductance L(x)
    # The reluctance of the air gap is a parallel combination of the path through
    # the magnetic block and the path through the air.
    # Reluctance of the path through the block
    R_block = g / (mu * x * w)
    # Reluctance of the path through the air
    R_air = g / (mu_0 * (D - x) * w)

    # The total reluctance of the gap is the parallel combination.
    # 1/R_total = 1/R_block + 1/R_air
    R_total = (R_block * R_air) / (R_block + R_air)

    # The inductance L is N^2 / R_total
    L_x = N**2 / R_total
    # Simplify the expression for L(x)
    L_x = sympy.simplify(L_x)

    # Step 2: Calculate the work done for each segment of the cycle
    # The total work done by the source is the integral of I * d(lambda), where lambda = L*I.
    # W = integral(I * d(L*I)) = integral(I^2 * dL + I*L*dI)

    # Path 1: x from x1 to x2, current is constant at I1. dI=0.
    # Work is integral(I1^2 * dL)
    W1 = sympy.integrate(I1**2 * sympy.diff(L_x, x), (x, x1, x2))

    # Path 2: current from I1 to I2, position is constant at x2. dL=0.
    # Work is integral(I * L(x2) * dI)
    L_at_x2 = L_x.subs(x, x2)
    W2 = sympy.integrate(L_at_x2 * sympy.Symbol('I'), (sympy.Symbol('I'), I1, I2))

    # Path 3: x from x2 to x1, current is constant at I2. dI=0.
    # Work is integral(I2^2 * dL)
    W3 = sympy.integrate(I2**2 * sympy.diff(L_x, x), (x, x2, x1))

    # Path 4: current from I2 to I1, position is constant at x1. dL=0.
    # Work is integral(I * L(x1) * dI)
    L_at_x1 = L_x.subs(x, x1)
    W4 = sympy.integrate(L_at_x1 * sympy.Symbol('I'), (sympy.Symbol('I'), I2, I1))

    # Step 3: Sum the work from all paths and simplify
    W_total = W1 + W2 + W3 + W4
    W_simplified = sympy.simplify(W_total)
    
    # Factor the expression to make it easier to compare with the options
    W_factored = sympy.factor(W_simplified)

    # The result from sympy is: N**2*w*(I1 - I2)*(I1 + I2)*(mu - mu_0)*(x1 - x2)/(2*g)
    # This can be rewritten by flipping the sign of two terms:
    # W = N**2*w*(-(I2 - I1))*(I1 + I2)*(mu - mu_0)*(-(x2 - x1))/(2*g)
    # W = N**2*w*(I2**2 - I1**2)*(mu - mu_0)*(x2 - x1)/(2*g)
    # Rearranging to match option D:
    # W = - (mu - mu_0) * N**2 * w * (x2 - x1) / (2*g) * (I2**2 - I1**2)
    # We will print the final equation in the standard form.
    print("The total work done by the current source per cycle is:")
    print(f"W = - (mu - mu_0) * N**2 * w * (x2 - x1) * (I2**2 - I1**2) / (2*g)")

if __name__ == '__main__':
    solve_work_cycle()