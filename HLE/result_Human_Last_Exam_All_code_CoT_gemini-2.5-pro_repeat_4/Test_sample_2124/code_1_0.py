def calculate_counter_term_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The calculation is based on the relationships between the counter-terms
    in the MS-bar scheme.
    """

    # Let C be the common factor g^2 / (32 * pi^2 * epsilon)
    # The coefficients of the counter-terms are defined relative to C.

    # 1. Define the coefficients for the primary counter-terms based on standard 1-loop results.
    # Coefficient for the fermion wave-function counter-term (delta Z_x)
    dZx_coeff = -1
    # Coefficient for the fermion mass counter-term (delta Z_mx)
    dZmx_coeff = -2
    # Coefficient for the coupling constant counter-term (delta_g)
    dg_coeff = -1

    print("Step 1: Define the coefficients of the counter-terms relative to a common factor C.")
    print(f"Coefficient for delta Z_x = {dZx_coeff}")
    print(f"Coefficient for delta Z_mx = {dZmx_coeff}")
    print(f"Coefficient for the coupling constant CT (delta_g) = {dg_coeff}\n")

    # 2. Calculate the coefficient for the full vertex counter-term (delta Z_g).
    # The problem states that the scalar field renormalization delta Z_phi = 0.
    # The relation is: delta Z_g = delta_g + delta Z_x + (1/2)*delta Z_phi
    # With delta Z_phi = 0, this simplifies to: delta Z_g = delta_g + delta Z_x
    dZg_coeff = dg_coeff + dZx_coeff

    print("Step 2: Calculate the coefficient for the full vertex counter-term delta Z_g.")
    print("The relation is delta Z_g = delta_g + delta Z_x (since delta Z_phi = 0).")
    print(f"Coefficient for delta Z_g = {dg_coeff} + ({dZx_coeff}) = {dZg_coeff}\n")

    # 3. Calculate the ratio R = delta Z_x / (delta Z_g + delta Z_mx)
    numerator = dZx_coeff
    denominator = dZg_coeff + dZmx_coeff
    
    # Avoid division by zero, though not expected here.
    if denominator == 0:
        result = "undefined"
    else:
        result = numerator / denominator

    print("Step 3: Calculate the ratio R = delta Z_x / (delta Z_g + delta Z_mx).")
    print("Plugging in the coefficients:")
    print(f"R = ({numerator}) / (({dZg_coeff}) + ({dZmx_coeff}))")
    print(f"R = ({numerator}) / ({denominator})")
    print(f"R = {result}\n")
    
    print("Final Answer:")
    print(f"The calculated ratio R is {result}")

calculate_counter_term_ratio()
<<<0.25>>>