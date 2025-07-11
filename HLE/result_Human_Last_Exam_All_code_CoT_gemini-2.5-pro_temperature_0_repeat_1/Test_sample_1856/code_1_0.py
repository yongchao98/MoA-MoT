import math

def calculate_decay_angle():
    """
    Calculates the angle of particle C in the lab frame.
    """
    # Given velocity of particle A
    beta_A = 0.95

    # Step 1: Calculate the Lorentz factor gamma for particle A
    gamma_A = 1 / math.sqrt(1 - beta_A**2)

    # The angle theta_C in the lab frame is given by the formula:
    # tan(theta_C) = p_transverse / p_longitudinal
    # p_transverse = p*_Cx = p*_C / sqrt(2)
    # p_longitudinal = gamma_A * (p*_Cz + beta_A * E*_C)
    # Using p*_Cz = p*_C / sqrt(2) and the ultra-relativistic approximation E*_C ≈ p*_C,
    # the p*_C terms cancel out, leaving:
    # tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))

    # Step 2: Calculate the terms in the formula
    val_1_over_sqrt_2 = 1 / math.sqrt(2)
    
    # Step 3: Calculate tan(theta_C)
    denominator_val = gamma_A * (val_1_over_sqrt_2 + beta_A)
    tan_theta_C = val_1_over_sqrt_2 / denominator_val

    # Step 4: Calculate the angle in degrees
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)

    # Print the breakdown of the calculation
    print("The formula for the tangent of the angle theta_C in the lab frame is:")
    print("tan(theta_C) = (1/sqrt(2)) / (gamma_A * (1/sqrt(2) + beta_A))")
    print("\nSubstituting the numerical values:")
    print(f"beta_A = {beta_A:.3f}")
    print(f"gamma_A = 1 / sqrt(1 - {beta_A:.3f}^2) = {gamma_A:.3f}")
    print(f"1/sqrt(2) = {val_1_over_sqrt_2:.3f}")
    
    print("\nPlugging these into the formula:")
    print(f"tan(theta_C) = {val_1_over_sqrt_2:.3f} / ({gamma_A:.3f} * ({val_1_over_sqrt_2:.3f} + {beta_A:.3f}))")
    print(f"tan(theta_C) = {val_1_over_sqrt_2:.3f} / {denominator_val:.3f}")
    print(f"tan(theta_C) = {tan_theta_C:.3f}")

    print("\nFinally, the angle in degrees is:")
    print(f"theta_C = arctan({tan_theta_C:.3f})")
    print(f"The angle is {theta_C_deg:.3f} degrees.")

if __name__ == '__main__':
    calculate_decay_angle()