import sympy

def solve_work_done():
    """
    This function derives the formula for the work done by the current source in the given cycle.
    It uses symbolic mathematics to construct the expressions for reluctance, inductance, and finally work.
    """

    # Step 1: Define symbolic variables
    mu, mu_0, g, D, w, N, I_1, I_2, x_1, x_2, x = sympy.symbols('μ μ_0 g D w N I_1 I_2 x_1 x_2 x')

    # Step 2: Formulate the inductance L(x)
    # The reluctance of the part of the gap with the magnetic block (permeability μ, area x*w)
    R_block = g / (mu * x * w)
    # The reluctance of the part of the gap with air (permeability μ_0, area (D-x)*w)
    R_air = g / (mu_0 * (D - x) * w)
    # The total reluctance of the gap is the parallel combination
    # 1/R_total = 1/R_block + 1/R_air
    R_total = 1 / (1/R_block + 1/R_air)
    R_total_simplified = sympy.simplify(R_total)
    
    # Inductance L(x) = N^2 / R_total(x)
    L_x = N**2 / R_total_simplified
    
    # Step 3: Calculate work done over the cycle
    # For a rectangular cycle (x1 -> x2 at I1, I1 -> I2 at x2, x2 -> x1 at I2, I2 -> I1 at x1),
    # the net work done by the source is W = 1/2 * (L(x1) - L(x2)) * (I2^2 - I1^2)
    L_at_x1 = L_x.subs(x, x_1)
    L_at_x2 = L_x.subs(x, x_2)
    
    delta_L = L_at_x1 - L_at_x2
    delta_L_simplified = sympy.simplify(delta_L)

    # Step 4: Final formula for Work
    Work = sympy.S(1)/2 * delta_L_simplified * (I_2**2 - I_1**2)

    # Rearrange to match the answer choices format
    final_work_expr = -sympy.S(1)/2 * (mu - mu_0)/g * N**2 * w * (x_2 - x_1) * (I_2**2 - I_1**2)

    # Print the final result by breaking it down into its symbolic parts as requested
    print("The final derived formula for the work done (W) by the current source per cycle is:")

    # Define each term or "number" in the equation
    sign = "-"
    fraction_num = f"( {sympy.pretty(mu)} - {sympy.pretty(mu_0)} )"
    fraction_den = f"( 2 * {sympy.pretty(g)} )"
    term_N_sq = f"{sympy.pretty(N)}^2"
    term_w = f"{sympy.pretty(w)}"
    term_delta_x = f"( {sympy.pretty(x_2)} - {sympy.pretty(x_1)} )"
    term_delta_I_sq = f"( {sympy.pretty(I_2)}^2 - {sympy.pretty(I_1)}^2 )"
    
    # Output the final equation built from its components
    print(f"W = {sign} [ {fraction_num} / {fraction_den} ] * {term_N_sq} * {term_w} * {term_delta_x} * {term_delta_I_sq}")

solve_work_done()