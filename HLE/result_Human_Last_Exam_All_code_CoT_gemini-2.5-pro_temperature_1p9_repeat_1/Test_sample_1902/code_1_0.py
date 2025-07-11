import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on the provided parameters.
    """
    # --- Input parameters from the problem statement ---
    # Photon beam size at focus (horizontal by vertical)
    beam_width_focus_mm = 0.3  # mm
    beam_height_focus_mm = 6.0   # mm
    # Ionization chamber current
    I_pA = 2.0                 # pA
    # Density of air in the ionization chamber
    rho_air_mg_cm3 = 1.293    # mg/cm^3
    # Length of ionization chamber
    L_chamber_cm = 15.1        # cm
    # Ratio of beam’s width at focus to the scan step of subject per exposure time
    t_total_s = 0.02           # s

    # --- Physical constants ---
    # Mean energy to create an ion pair in air / elementary charge (W_air/e)
    W_over_e_J_C = 33.97  # J/C

    # --- Step 1: Convert units to a consistent system (SI: meters, kg, seconds, Coulombs) ---
    print("Converting all inputs to standard SI units (m, kg, s, A)...\n")
    beam_width_focus_m = beam_width_focus_mm / 1000.0
    beam_height_focus_m = beam_height_focus_mm / 1000.0
    L_chamber_m = L_chamber_cm / 100.0
    I_A = I_pA * 1e-12  # Convert picoamperes (pA) to amperes (A or C/s)
    
    # Convert density from mg/cm^3 to kg/m^3
    # 1 mg/cm^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
    rho_air_kg_m3 = rho_air_mg_cm3

    # --- Step 2: Calculate the irradiated volume of the ionization chamber ---
    V_irradiated_m3 = beam_width_focus_m * beam_height_focus_m * L_chamber_m

    # --- Step 3: Calculate the mass of the air inside the irradiated volume ---
    m_air_kg = V_irradiated_m3 * rho_air_kg_m3

    # --- Step 4: Calculate the absorbed dose rate in air (D_dot_air) ---
    # Dose Rate (Gray/s) = Energy deposited per second (J/s) / mass (kg)
    # Energy deposited per second (J/s) = Current (C/s) * W/e (J/C)
    D_dot_air_Gy_s = (I_A * W_over_e_J_C) / m_air_kg
    
    # Per the problem, assume dose rate in tissue is equal to dose rate in air
    D_dot_tissue_Gy_s = D_dot_air_Gy_s

    # --- Step 5: Calculate the cumulative surface dose ---
    # Cumulative Dose (Gy) = Dose Rate (Gy/s) * Total Exposure Time (s)
    D_cumulative_Gy = D_dot_tissue_Gy_s * t_total_s

    # --- Step 6: Convert results to milligray (mGy) for readability ---
    D_cumulative_mGy = D_cumulative_Gy * 1000.0
    D_dot_tissue_mGy_s = D_dot_tissue_Gy_s * 1000.0

    # --- Print the breakdown and final calculation ---
    print("--- Calculation Breakdown ---")
    print(f"1. Irradiated Volume: {beam_width_focus_m:.3e} m * {beam_height_focus_m:.3e} m * {L_chamber_m:.3e} m = {V_irradiated_m3:.4e} m^3")
    print(f"2. Mass of Irradiated Air: {V_irradiated_m3:.4e} m^3 * {rho_air_kg_m3:.3f} kg/m^3 = {m_air_kg:.4e} kg")
    print(f"3. Tissue Dose Rate: ({I_A:.2e} C/s * {W_over_e_J_C:.2f} J/C) / {m_air_kg:.4e} kg = {D_dot_tissue_mGy_s:.6f} mGy/s")
    print(f"4. Total Exposure Time for a single point: {t_total_s:.2f} s")
    print("\n--- Final Cumulative Dose Calculation ---")
    print("Cumulative Dose (mGy) = Dose Rate (mGy/s) * Exposure Time (s)")
    print(f"{D_cumulative_mGy:.6f} mGy = {D_dot_tissue_mGy_s:.6f} mGy/s * {t_total_s:.2f} s")
    
    # --- Final Answer in specified format ---
    print(f"\n<<<{D_cumulative_mGy:.6f}>>>")

if __name__ == '__main__':
    calculate_surface_dose()