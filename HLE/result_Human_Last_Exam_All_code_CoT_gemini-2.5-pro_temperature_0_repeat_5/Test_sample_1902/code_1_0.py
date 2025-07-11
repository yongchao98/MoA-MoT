import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- Given Information ---
    # Beam dimensions at focus
    beam_width_focus_mm = 0.3  # mm
    beam_height_focus_mm = 6.0   # mm
    # Ionization chamber properties
    chamber_length_cm = 15.1 # cm
    current_pA = 2.0         # pA
    air_density_mg_cm3 = 1.293 # mg/cm^3
    # Exposure time for a point on the subject
    exposure_time_s = 0.02   # s
    # Physical constant: mean energy to create an ion pair in air per unit charge
    W_air_over_e_J_C = 33.97 # J/C

    # --- Unit Conversions to SI units (meters, kg, seconds, Amperes) ---
    beam_width_focus_m = beam_width_focus_mm / 1000.0
    beam_height_focus_m = beam_height_focus_mm / 1000.0
    chamber_length_m = chamber_length_cm / 100.0
    current_A = current_pA * 1e-12
    # Convert density from mg/cm^3 to kg/m^3
    # 1 mg/cm^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
    air_density_kg_m3 = air_density_mg_cm3

    # --- Step 1: Calculate the irradiated volume of air ---
    beam_area_m2 = beam_width_focus_m * beam_height_focus_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m

    # --- Step 2: Calculate the mass of the irradiated air ---
    irradiated_mass_kg = irradiated_volume_m3 * air_density_kg_m3

    # --- Step 3: Calculate the dose rate in Gray per second (J/kg/s) ---
    # Dose Rate = (Current / Mass) * (Energy per charge)
    dose_rate_Gy_s = (current_A / irradiated_mass_kg) * W_air_over_e_J_C

    # --- Step 4: Calculate the cumulative dose in Gray (Gy) ---
    # Cumulative Dose = Dose Rate * Exposure Time
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

    # --- Convert final dose to milliGray (mGy) for readability ---
    cumulative_dose_mGy = cumulative_dose_Gy * 1000

    # --- Print the calculation steps with values ---
    print("Calculating the cumulative surface dose using the formula:")
    print("Dose = (Current / (Beam_Area * Chamber_Length * Air_Density)) * (W_air/e) * Exposure_Time\n")

    print("--- Input Values (in SI units) ---")
    print(f"Current (I): {current_A:.2e} A")
    print(f"Beam Area (A): {beam_area_m2:.2e} m^2")
    print(f"Chamber Length (L): {chamber_length_m:.4f} m")
    print(f"Air Density (ρ): {air_density_kg_m3:.4f} kg/m^3")
    print(f"Energy per Charge (W_air/e): {W_air_over_e_J_C:.2f} J/C")
    print(f"Exposure Time (t): {exposure_time_s:.2f} s\n")

    print("--- Calculation Breakdown ---")
    print(f"Irradiated Mass (m = A * L * ρ): {irradiated_mass_kg:.4e} kg")
    print(f"Dose Rate (D_rate = (I / m) * (W_air/e)): {dose_rate_Gy_s:.4e} Gy/s")
    print(f"Cumulative Dose (Dose = D_rate * t): {cumulative_dose_Gy:.4e} Gy\n")

    print("--- Final Answer ---")
    print(f"The cumulative surface dose is {cumulative_dose_mGy:.5f} mGy.")

if __name__ == '__main__':
    calculate_surface_dose()
<<<0.00387>>>