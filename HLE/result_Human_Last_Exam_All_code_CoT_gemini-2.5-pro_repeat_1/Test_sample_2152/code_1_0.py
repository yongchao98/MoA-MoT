import numpy as np

def calculate_cross_section_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering
    by a magnetic monopole and a magnetic dipole.
    """
    # --- Constants and Given Values ---
    # Particle properties
    v_c_ratio = 1/100  # v = c/100

    # Monopole properties
    e_m_e_ratio = 1/16  # e_m = e/16

    # Dipole properties
    mu_muB_ratio = 25  # mu = 25 * mu_B

    # Scattering angle
    theta = np.pi / 30  # radians

    # Physical constants
    c = 299792458  # Speed of light in m/s

    # Particle speed
    v = c * v_c_ratio

    # --- Derivation of the Ratio ---
    # The ratio of the differential cross-sections is:
    # R = ( (mu_0 * e * e_m) / (8 * pi * p * sin^2(theta/2)) )^2 / ( (mu_0 * e * mu * cot(theta/2)) / (4 * pi * h_bar * sqrt(3)) )^2
    # After simplification, this becomes:
    # R = ( (e_m * h_bar * sqrt(3)) / (2 * p * mu * sin(theta/2) * cos(theta/2)) )^2
    # R = ( (e_m * h_bar * sqrt(3)) / (p * mu * sin(theta)) )^2
    
    # Substitute e_m = (e/16), mu = 25 * (e*h_bar / (2*m_e)), and p = m_e * v
    # The terms e, h_bar, and m_e cancel out, leading to:
    # R = ( ((1/16) * sqrt(3)) / (v * (25/2) * sin(theta)) )^2
    # R = ( (2 * sqrt(3)) / (16 * 25 * v * sin(theta)) )^2
    # R = ( sqrt(3) / (200 * v * sin(theta)) )^2
    # Substitute v = c/100:
    # R = ( sqrt(3) / (200 * (c/100) * sin(theta)) )^2
    # R = ( sqrt(3) / (2 * c * sin(theta)) )^2

    # --- Final Calculation ---
    numerator = np.sqrt(3)
    denominator = 2 * c * np.sin(theta)
    
    ratio = (numerator / denominator)**2

    # --- Output Results ---
    print("This script calculates the ratio of differential cross-sections for a charged particle scattered by a magnetic monopole versus a magnetic dipole.")
    print("\nGiven values:")
    print(f"  Particle speed, v = c/100 = {v:.2e} m/s")
    print(f"  Monopole charge, e_m = e/16")
    print(f"  Dipole moment, mu = 25 * mu_B")
    print(f"  Scattering angle, theta = pi/30 = {np.degrees(theta):.2f} degrees")
    
    print("\nThe formula for the ratio, after simplification, is:")
    print("  Ratio = (sqrt(3) / (2 * c * sin(theta)))^2")

    print("\nPlugging in the numbers:")
    final_eq = f"  Ratio = (sqrt(3) / (2 * {c} * sin(pi/30)))^2"
    print(final_eq)
    final_eq_eval = f"  Ratio = ({numerator:.6f} / (2 * {c} * {np.sin(theta):.6f}))^2"
    print(final_eq_eval)
    final_eq_result = f"  Ratio = ({numerator / denominator:.6e})^2"
    print(final_eq_result)

    print(f"\nThe final calculated ratio is: {ratio:.4e}")
    
    # Returning the value for the final answer block
    return ratio

final_answer = calculate_cross_section_ratio()
# The final answer is requested in a specific format.
# The print statements above provide the detailed explanation.
# The following line will not be printed in the final output but is used to generate the answer block.
# print(f'<<<{final_answer:.4e}>>>')