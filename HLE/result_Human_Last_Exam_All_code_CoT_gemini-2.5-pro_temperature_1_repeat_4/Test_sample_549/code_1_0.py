import math

def calculate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity for an electron in a 3D
    bulk semiconductor due to the weak localization effect.
    """
    # 1. Define constants and typical material parameters
    # Fundamental Constants
    e = 1.60217663e-19  # Elementary charge in Coulombs
    hbar = 1.05457182e-34 # Reduced Planck constant in J·s
    m_e = 9.1093837e-31   # Electron rest mass in kg
    
    # Typical parameters for a doped semiconductor at low temperature
    n = 1e23              # Electron concentration in m^-3 (e.g., 10^17 cm^-3)
    m_star = 0.1 * m_e    # Effective mass (e.g., 0.1 * m_e for some III-V semiconductors)
    tau = 1e-13           # Elastic mean free time (momentum relaxation time) in seconds
    tau_phi = 1e-11       # Phase coherence time in seconds (typically >> tau)

    print("--- Input Parameters ---")
    print(f"Electron concentration (n): {n:.1e} m^-3")
    print(f"Effective mass (m*): {m_star:.2e} kg")
    print(f"Elastic scattering time (τ): {tau:.1e} s")
    print(f"Phase coherence time (τφ): {tau_phi:.1e} s\n")

    # 2. Calculate Fermi Energy (E_F)
    E_F = (hbar**2 / (2 * m_star)) * (3 * math.pi**2 * n)**(2/3)
    print("--- Step 1: Calculate Fermi Energy (E_F) ---")
    print(f"E_F = (ħ² / (2 * m*)) * (3 * π² * n)^(2/3)")
    print(f"E_F = {E_F:.3e} J (or {E_F / e * 1000:.2f} meV)\n")

    # 3. Calculate Fermi Velocity (v_F)
    v_F = math.sqrt(2 * E_F / m_star)
    print("--- Step 2: Calculate Fermi Velocity (v_F) ---")
    print(f"v_F = sqrt(2 * E_F / m*)")
    print(f"v_F = {v_F:.3e} m/s\n")

    # 4. Calculate Diffusion Constant (D)
    D = (1/3) * v_F**2 * tau
    print("--- Step 3: Calculate Diffusion Constant (D) ---")
    print(f"D = (1/3) * v_F² * τ")
    print(f"D = {D:.3e} m²/s\n")

    # 5. Calculate Phase Coherence Length (L_phi)
    L_phi = math.sqrt(D * tau_phi)
    print("--- Step 4: Calculate Phase Coherence Length (Lφ) ---")
    print(f"Lφ = sqrt(D * τφ)")
    print(f"Lφ = {L_phi:.3e} m (or {L_phi * 1e9:.2f} nm)\n")

    # 6. Calculate the Quantum Correction to Conductivity (Δσ)
    delta_sigma = -(e**2) / (2 * math.pi**2 * hbar * L_phi)
    print("--- Step 5: Final Calculation for Quantum Correction (Δσ) ---")
    print("The formula for the quantum correction to conductivity in 3D is:")
    print("Δσ = -e² / (2 * π² * ħ * Lφ)\n")

    print("Substituting the calculated and constant values:")
    print(f"Δσ = -({e:.4e} C)² / (2 * π² * {hbar:.4e} J·s * {L_phi:.4e} m)")
    
    # Final value
    prefactor_val = (e**2) / (2 * math.pi**2 * hbar)
    final_val = -prefactor_val / L_phi
    
    print(f"Δσ = -({prefactor_val:.4e} S·m) / ({L_phi:.4e} m)")
    print(f"Δσ = {final_val:.2f} S/m")

    # Optional: Compare with classical Drude conductivity
    sigma_0 = (n * e**2 * tau) / m_star
    print(f"\nFor context, the classical Drude conductivity (σ₀) is {sigma_0:.2f} S/m.")
    print(f"The quantum correction is {final_val / sigma_0 * 100:.2f}% of the classical value.")
    
    return final_val

if __name__ == '__main__':
    final_answer = calculate_conductivity_correction()
    # The final answer is wrapped according to the required format.
    # print(f"\n<<< {final_answer:.2f} >>>")
    # Let's directly print the value
    # print(f"<<<{final_answer:.1f}>>>")

calculate_conductivity_correction()
<<< -128.49 >>>