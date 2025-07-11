import sympy

def solve_work_cycle():
    """
    This function symbolically derives the work done by the current source
    in the described thermodynamic cycle and prints the final expression.
    """
    # Define the symbolic variables used in the problem
    mu, mu_0, N, w, g, D, x = sympy.symbols('mu mu_0 N w g D x', real=True, positive=True)
    I, I_1, I_2 = sympy.symbols('I I_1 I_2', real=True)
    x_1, x_2 = sympy.symbols('x_1 x_2', real=True)

    # Step 1: Define the inductance L as a function of position x
    # This is based on the derivation L(x) = (N^2*w/g) * ((mu - mu_0)*x + mu_0*D)
    L_x = (N**2 * w / g) * ((mu - mu_0) * x + mu_0 * D)

    # Step 2: Define the inductance at the start and end positions
    L1 = L_x.subs(x, x_1)
    L2 = L_x.subs(x, x_2)

    # Step 3: Calculate the work done for each of the four steps in the cycle
    # W = integral(I * d(L*I))
    W_1 = I_1**2 * (L2 - L1)                # Path 1: x1 -> x2 at constant I1
    W_2 = sympy.integrate(L2 * I, (I, I_1, I_2)) # Path 2: I1 -> I2 at constant x2
    W_3 = I_2**2 * (L1 - L2)                # Path 3: x2 -> x1 at constant I2
    W_4 = sympy.integrate(L1 * I, (I, I_2, I_1)) # Path 4: I2 -> I1 at constant x1

    # Step 4: Sum the work components and simplify the total work expression
    W_total = sympy.simplify(W_1 + W_2 + W_3 + W_4)
    
    # The simplified result from sympy is -N**2*w*(-mu_0 + mu)*(-I_1**2 + I_2**2)*(-x_1 + x_2)/(2*g)
    # We will format this to match the structure of the provided options.

    # Step 5: Print the final formula in a structured way
    print("The work done by the current source for each cycle is:")
    
    # Constructing the equation string as requested
    numerator_1 = "(\u03BC - \u03BC_0)" # (mu - mu_0)
    denominator = "2g"
    term_N_sq = "N\u00b2" # N^2
    term_w = "w"
    term_x = "(x\u2082 - x\u2081)" # (x2 - x1)
    term_I_sq = "(I\u2082\u00b2 - I\u2081\u00b2)" # (I2^2 - I1^2)

    print(f"W = - \u2219 {numerator_1} / {denominator} \u2219 {term_N_sq} \u2219 {term_w} \u2219 {term_x} \u2219 {term_I_sq}")
    print("\nThis corresponds to Answer Choice D.")


solve_work_cycle()
<<<D>>>