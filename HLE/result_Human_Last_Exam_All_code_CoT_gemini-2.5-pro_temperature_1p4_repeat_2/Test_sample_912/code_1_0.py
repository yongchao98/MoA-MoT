import sympy

def solve_work_cycle():
    """
    This function symbolically derives the work done by the current source
    over a rectangular cycle in the I-x plane for a given electromagnetic actuator.
    """
    # Define the symbolic variables
    mu_0, mu, N, g, D, w, x = sympy.symbols('mu_0 mu N g D w x', real=True, positive=True)
    x1, x2 = sympy.symbols('x1 x2', real=True, positive=True)
    I, I1, I2 = sympy.symbols('I I1 I2', real=True)

    # Step 1: Define the inductance L(x)
    # Reluctance of the part of the gap filled by the block
    R_block = g / (mu * x * w)
    # Reluctance of the part of the gap filled by air
    R_air = g / (mu_0 * (D - x) * w)
    
    # Total reluctance of the gap (parallel paths)
    R_gap = 1 / (1/R_block + 1/R_air)

    # Inductance L(x) = N^2 / R_total (where R_total = R_gap)
    L = N**2 / R_gap
    
    # Simplify the expression for inductance
    L_simplified = sympy.simplify(L)
    
    # Step 2: Calculate work for each step of the cycle
    L_x1 = L_simplified.subs(x, x1)
    L_x2 = L_simplified.subs(x, x2)
    
    # Step A: x from x1 to x2, current is I1
    W_A = I1**2 * (L_x2 - L_x1)
    
    # Step B: current from I1 to I2, position is x2
    W_B = sympy.integrate(L_x2 * I, (I, I1, I2))
    
    # Step C: x from x2 to x1, current is I2
    W_C = I2**2 * (L_x1 - L_x2)
    
    # Step D: current from I2 to I1, position is x1
    W_D = sympy.integrate(L_x1 * I, (I, I2, I1))

    # Step 3: Sum the work from all four steps
    W_total = W_A + W_B + W_C + W_D
    
    # Step 4: Simplify the total work expression
    W_final = sympy.factor(W_total)

    # Print the final result
    print("The work done by the current source for each cycle is:")
    # The pretty print function displays the equation in a more readable format.
    # The equation shows each symbolic "number" or variable.
    sympy.pprint(sympy.Eq(sympy.Symbol('W'), W_final), use_unicode=False)

if __name__ == '__main__':
    solve_work_cycle()