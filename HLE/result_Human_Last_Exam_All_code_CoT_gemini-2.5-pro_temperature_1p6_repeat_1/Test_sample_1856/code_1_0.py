import math

def calculate_decay_angle():
    """
    Calculates the angle between particle C and its parent particle A in the lab frame.
    """
    # Step 1: Define the given constant.
    beta_A = 0.95

    # Step 2: Calculate the Lorentz factor for particle A.
    # The formula is gamma = 1 / sqrt(1 - beta^2).
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # Step 3: Derive the formula for the angle in the lab frame.
    # The transformation of momentum from A's rest frame (starred) to the lab frame is:
    # p_Cx = p_Cx* = p_C* / sqrt(2)
    # p_Cz = gamma_A * (p_Cz* + beta_A * E_C*)
    #
    # Using the ultra-relativistic approximation E_C* = |p_C*| = p_C*, we get:
    # p_Cz = gamma_A * (p_C*/sqrt(2) + beta_A * p_C*)
    #
    # The angle theta_C is given by tan(theta_C) = p_Cx / p_Cz.
    # The unknown momentum p_C* cancels out, leading to:
    # tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))
    # which simplifies to tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))

    # Step 4: Perform the numerical calculation.
    sqrt2 = math.sqrt(2)
    denominator_val = gamma_A * (1 + sqrt2 * beta_A)
    tan_theta_C = 1 / denominator_val
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)

    # Step 5: Print the detailed explanation and the result.
    print("The angle of particle C in the lab frame is found by a Lorentz transformation.")
    print("The formula for the tangent of the angle, tan(theta_C), simplifies to:")
    print("tan(theta_C) = 1 / (gamma_A * (1 + sqrt(2) * beta_A))")
    print("\n--- Calculation Steps ---")
    print(f"Given beta_A = {beta_A}")
    print(f"1. Calculate the Lorentz factor, gamma_A:")
    print(f"   gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.5f}")
    
    # Store intermediate calculation for printing
    term_in_parenthesis = 1 + sqrt2 * beta_A
    print("\n2. Calculate the denominator of the tangent formula:")
    print(f"   Denominator = gamma_A * (1 + sqrt(2) * beta_A)")
    print(f"   Denominator = {gamma_A:.5f} * (1 + {sqrt2:.5f} * {beta_A})")
    print(f"   Denominator = {gamma_A:.5f} * {term_in_parenthesis:.5f} = {denominator_val:.5f}")

    print("\n3. Calculate the tangent of the angle:")
    print(f"   tan(theta_C) = 1 / {denominator_val:.5f} = {tan_theta_C:.5f}")

    print("\n4. Calculate the angle in degrees:")
    print(f"   theta_C = arctan({tan_theta_C:.5f}) = {theta_C_deg:.5f} degrees")

    print("\n--- Final Answer ---")
    final_angle = round(theta_C_deg, 3)
    print("The final equation with the numerical values is:")
    print(f"Angle = arctan(1 / ({gamma_A:.3f} * (1 + {sqrt2:.3f} * {beta_A})))")
    print(f"The angle between particle C and its parent particle A is {final_angle} degrees.")
    
    # Final answer in the required format
    print(f"\n<<<{final_angle}>>>")

calculate_decay_angle()