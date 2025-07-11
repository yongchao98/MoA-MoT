import math

def calculate_lab_frame_angle():
    """
    Calculates the angle of particle C in the lab frame.

    This function implements the relativistic angle transformation to find the
    angle of a decay product in the lab frame, given its properties in the
    parent particle's rest frame.
    """

    # --- Step 1: Define given constants ---
    # Velocity of the parent particle A
    beta_A = 0.95
    
    # In the rest frame of A, the momentum components of C are given as
    # p_Cx* = P_C* / sqrt(2) and p_Cz* = P_C* / sqrt(2).
    # The angle theta_C* in the rest frame is given by tan(theta_C*) = p_Cx* / p_Cz*.
    # tan(theta_C*) = (P_C*/sqrt(2)) / (P_C*/sqrt(2)) = 1
    # This means theta_C* is 45 degrees.
    theta_C_star_deg = 45.0
    theta_C_star_rad = math.radians(theta_C_star_deg)

    # --- Step 2: Calculate the Lorentz factor gamma ---
    gamma_A = 1.0 / math.sqrt(1 - beta_A**2)

    # --- Step 3: Use the angle transformation formula ---
    # The formula is: tan(theta_C) = sin(theta_C*) / (gamma_A * (cos(theta_C*) + beta_A))
    # We will calculate each part of this equation.
    
    # Values of sin and cos for the angle in the rest frame
    sin_theta_star = math.sin(theta_C_star_rad)
    cos_theta_star = math.cos(theta_C_star_rad)
    
    # Denominator of the expression
    denominator = gamma_A * (cos_theta_star + beta_A)
    
    # Value for tan(theta_C)
    tan_theta_C = sin_theta_star / denominator

    # --- Step 4: Calculate the final angle in degrees ---
    theta_C_rad = math.atan(tan_theta_C)
    theta_C_deg = math.degrees(theta_C_rad)
    
    # --- Step 5: Print the breakdown of the calculation and the final answer ---
    print("Final Equation for the angle theta_C in the lab frame:")
    print("theta_C = arctan( sin(theta_C*) / (gamma_A * (cos(theta_C*) + beta_A)) )\n")

    print("--- Calculation Breakdown ---")
    print(f"Parent particle A velocity, beta_A              = {beta_A}")
    print(f"Particle C angle in rest frame, theta_C*        = {theta_C_star_deg} degrees")
    print(f"Lorentz factor for A, gamma_A                 = {gamma_A:.6f}")
    print(f"sin({theta_C_star_deg} deg)                                   = {sin_theta_star:.6f}")
    print(f"cos({theta_C_star_deg} deg)                                   = {cos_theta_star:.6f}")
    print(f"Denominator term: gamma_A * (cos_theta* + beta_A) = {denominator:.6f}")
    print(f"tan(theta_C) = sin_theta* / Denominator       = {tan_theta_C:.6f}\n")
    
    print("--- Final Result ---")
    print("The angle between particle C and particle A in the lab frame is:")
    print(f"theta_C = {theta_C_deg:.3f} degrees")

if __name__ == "__main__":
    calculate_lab_frame_angle()
<<<7.588>>>