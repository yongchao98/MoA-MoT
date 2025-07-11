import numpy as np
from scipy.constants import mu_0

def calculate_magnetic_field():
    """
    Calculates the magnetic field B_z outside an infinite stack of superconducting strips.
    """
    # --- Define physical parameters (in SI units) ---
    # Critical current density (A/m^2)
    Jc = 1.0e10
    # Strip thickness (m)
    d = 1.0e-6
    # Strip half-width (m)
    w = 100.0e-6
    # Stacking interval (m)
    D = 200.0e-6
    # Applied magnetic field (A/m)
    Ha = 5000.0
    # Position to calculate the field (x, z) in meters
    x = 300.0e-6
    z = 50.0e-6

    # --- Step 1: Calculate characteristic field H0 ---
    H0 = Jc * d / np.pi

    # Check if the condition Ha > H0 is met
    if Ha <= H0:
        print("Warning: The condition Ha > H0 is not met.")
        print(f"Ha = {Ha:.2f} A/m, H0 = {H0:.2f} A/m")
        # The derivation is not valid in this regime.
        # We can still proceed to show the calculation, but the result is not physically accurate.

    # --- Step 2: Calculate the magnitude of the magnetic moment per unit length (m_0) ---
    # m_0 = pi * w^2 * H0 * tanh^2(Ha / H0)
    tanh_term = np.tanh(Ha / H0)**2
    m0_factor = np.pi * w**2 * H0
    m0 = m0_factor * tanh_term
    
    # --- Step 3: Calculate the geometric factor from the field of the dipole array ---
    # G = (sinh^2(pi*x/D) - sin^2(pi*z/D)) / (sinh^2(pi*x/D) + sin^2(pi*z/D))^2
    u = np.pi * x / D
    v = np.pi * z / D
    sinh_u_sq = np.sinh(u)**2
    sin_v_sq = np.sin(v)**2
    
    numerator = sinh_u_sq - sin_v_sq
    denominator = (sinh_u_sq + sin_v_sq)**2
    geometric_factor = numerator / denominator

    # --- Step 4: Calculate the induced magnetic field Hz_ind ---
    # Hz_ind = (pi * m0 / (2*D^2)) * G
    prefactor_H = np.pi / (2 * D**2)
    Hz_ind = prefactor_H * m0 * geometric_factor
    
    # --- Step 5: Calculate the total magnetic induction Bz ---
    # Bz = mu_0 * (Ha + Hz_ind)
    Hz_total = Ha + Hz_ind
    Bz_total = mu_0 * Hz_total
    
    # --- Print the results and the equation with numbers ---
    print("Derived Expression for the Magnetic Field B_z(x, z):")
    print("B_z = mu_0 * [ Ha + (pi^2 * w^2 * H0 / (2*D^2)) * tanh^2(Ha / H0) * G ]")
    print("where G = (sinh^2(pi*x/D) - sin^2(pi*z/D)) / (sinh^2(pi*x/D) + sin^2(pi*z/D))^2\n")
    
    print("--- Calculation with given values ---")
    print("Input Parameters:")
    print(f"  Ha = {Ha} A/m")
    print(f"  Jc = {Jc:.1e} A/m^2")
    print(f"  d = {d:.1e} m")
    print(f"  w = {w:.1e} m")
    print(f"  D = {D:.1e} m")
    print(f"  (x, z) = ({x:.1e} m, {z:.1e} m)\n")
    
    print("Intermediate Values:")
    print(f"  H0 = Jc*d/pi = {H0:.4f} A/m")
    print(f"  Ha/H0 = {Ha/H0:.4f}")
    print(f"  tanh^2(Ha/H0) = {tanh_term:.4f}")
    print(f"  Geometric Factor G = {geometric_factor:.4f}\n")

    print("Final Equation with Numbers:")
    # This shows the structure of the final calculation
    # B_z = mu_0 * ( Ha + Hz_ind )
    print(f"B_z = {mu_0:.2e} * ( {Ha:.2f} + {Hz_ind:.2f} )")
    print(f"B_z = {mu_0:.2e} * ( {Hz_total:.2f} )\n")

    print("--- Final Result ---")
    print(f"The magnetic field at (x={x*1e6:.0f} um, z={z*1e6:.0f} um) is B_z = {Bz_total:.6f} T")
    
    # Return the final numerical answer as requested by the format
    return Bz_total

# Execute the function and capture the final answer
final_answer = calculate_magnetic_field()
print(f"\n<<<{final_answer}>>>")
