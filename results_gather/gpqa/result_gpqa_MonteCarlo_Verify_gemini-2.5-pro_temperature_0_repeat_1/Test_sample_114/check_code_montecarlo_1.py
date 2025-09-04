import math

def check_electron_positron_annihilation():
    """
    This function verifies the answer to the electron-positron annihilation problem.
    It calculates the angle between the two resulting photons based on conservation
    of energy and momentum and compares it to the provided answer.
    """
    
    # --- Given values from the problem ---
    gamma_electron = 4.0
    gamma_positron = 2.0
    
    # The answer from option C to be checked
    answer_angle_deg = 138.0

    # --- Calculations based on physical principles (in natural units m_e=1, c=1) ---

    # 1. Conservation of Energy
    # Total initial energy E = γ_e * m_e + γ_p * m_e. With m_e=1, E_initial = γ_e + γ_p
    E_initial = gamma_electron + gamma_positron
    
    # The two photons have equal energy, so E_final = 2 * E_photon
    # By conservation, E_initial = E_final, so E_photon = E_initial / 2
    E_photon = E_initial / 2.0
    
    # 2. Conservation of Momentum
    # Relativistic momentum p = sqrt(γ^2 - 1)
    # Electron moves from left (+x direction), positron from right (-x direction)
    p_electron_x = math.sqrt(gamma_electron**2 - 1)
    p_positron_x = -math.sqrt(gamma_positron**2 - 1)
    
    # Total initial momentum is purely in the x-direction
    p_initial_x = p_electron_x + p_positron_x
    
    # For a photon, momentum magnitude p equals its energy E. So, p_photon = E_photon.
    p_photon = E_photon
    
    # The photons are emitted at angles +θ and -θ relative to the x-axis.
    # The final momentum in the y-direction is p_photon*sin(θ) - p_photon*sin(θ) = 0, which is conserved.
    # The final momentum in the x-direction is p_photon*cos(θ) + p_photon*cos(θ)
    # p_final_x = 2 * p_photon * cos(θ)
    
    # 3. Solve for the angle θ
    # By conservation of momentum, p_initial_x = p_final_x
    # p_initial_x = 2 * p_photon * cos(θ)
    # Therefore, cos(θ) = p_initial_x / (2 * p_photon)
    
    cos_theta = p_initial_x / (2.0 * p_photon)
    
    # Check if the value for arccos is valid
    if not -1.0 <= cos_theta <= 1.0:
        return f"Calculation Error: cos(θ) is {cos_theta:.4f}, which is outside the valid range [-1, 1]."

    # Calculate θ in radians
    theta_rad = math.acos(cos_theta)
    
    # The question asks for the angle *between* the two photons, which is 2 * θ
    total_angle_rad = 2.0 * theta_rad
    
    # Convert the final angle to degrees
    calculated_angle_deg = math.degrees(total_angle_rad)
    
    # 4. Verify the answer
    # Check if the calculated angle is close to the answer from option C.
    # A tolerance of 1 degree is sufficient given the options.
    if math.isclose(calculated_angle_deg, answer_angle_deg, abs_tol=1.0):
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {answer_angle_deg} degrees, but the "
                f"calculated angle is {calculated_angle_deg:.2f} degrees. The derivation "
                f"in the provided answer is correct, but the final choice should be the "
                f"closest integer, which is 138. The provided answer is correct, but this "
                f"checker code indicates a mismatch if the answer was different. "
                f"Calculated value: {calculated_angle_deg:.2f}, Answer value: {answer_angle_deg}")

# Run the check
result = check_electron_positron_annihilation()
print(result)