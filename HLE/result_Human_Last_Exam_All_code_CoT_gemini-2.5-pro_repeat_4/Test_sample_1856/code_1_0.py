import math

def calculate_decay_angle():
    """
    Calculates the angle of a decay product in the lab frame based on its
    emission in the parent particle's rest frame.
    """
    # Given velocity of particle A in the lab frame
    beta_A = 0.95

    # 1. Calculate the Lorentz factor gamma_A
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # 2. Determine the angle theta* in the rest frame of A.
    # The four-momentum C*_mu = (E*_C, P*_C/sqrt(2), 0, P*_C/sqrt(2)) gives:
    # p*_z = P*_C / sqrt(2) and p*_x = P*_C / sqrt(2).
    # The total momentum P*_C = sqrt(p*_x^2 + p*_y^2 + p*_z^2) = P*_C.
    # cos(theta*) = p*_z / P*_C = 1/sqrt(2)
    # sin(theta*) = p*_x / P*_C = 1/sqrt(2)
    cos_theta_star = 1 / math.sqrt(2)
    sin_theta_star = 1 / math.sqrt(2)
    theta_star_deg = math.degrees(math.acos(cos_theta_star))

    # 3. Use the relativistic aberration formula to find the angle theta_C in the lab frame.
    # tan(theta_C) = sin(theta*) / (gamma_A * (cos(theta*) + beta_A))
    
    # 4. Calculate the terms of the equation numerically.
    tan_theta_C_numerator = sin_theta_star
    denominator_in_parentheses = cos_theta_star + beta_A
    tan_theta_C_denominator = gamma_A * denominator_in_parentheses
    tan_theta_C = tan_theta_C_numerator / tan_theta_C_denominator

    # 5. Calculate the final angle in degrees.
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)

    # Print the step-by-step calculation as requested.
    print("The goal is to find the angle of particle C with respect to the beam axis (z-axis) in the lab frame.")
    print("\nStep 1: Define the kinematics.")
    print(f"The parent particle A moves with velocity beta_A = {beta_A}.")
    print(f"The corresponding Lorentz factor is gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.4f}.")
    print(f"In A's rest frame, particle C is emitted at an angle of theta* = {theta_star_deg:.1f} degrees.")
    
    print("\nStep 2: Use the angle transformation formula.")
    print("The lab frame angle (theta_C) is found using: tan(theta_C) = sin(theta*) / (gamma_A * (cos(theta*) + beta_A))")
    
    print("\nStep 3: Plug in the numbers into the equation.")
    print(f"tan(theta_C) = {sin_theta_star:.4f} / ({gamma_A:.4f} * ({cos_theta_star:.4f} + {beta_A}))")
    print(f"tan(theta_C) = {sin_theta_star:.4f} / ({gamma_A:.4f} * {denominator_in_parentheses:.4f})")
    print(f"tan(theta_C) = {sin_theta_star:.4f} / {tan_theta_C_denominator:.4f}")
    print(f"tan(theta_C) = {tan_theta_C:.4f}")

    print("\nStep 4: Solve for the angle.")
    print(f"theta_C = arctan({tan_theta_C:.4f})")
    print(f"The final angle in the lab frame is: {theta_C_deg:.3f} degrees.")

calculate_decay_angle()