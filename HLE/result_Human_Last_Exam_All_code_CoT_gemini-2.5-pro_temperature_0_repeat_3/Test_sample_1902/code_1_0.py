import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- Step 1: Define constants and given parameters ---

    # Given parameters from the problem description
    beam_width_h_mm = 0.3  # Focused beam horizontal width in mm
    beam_width_v_mm = 6.0  # Focused beam vertical width in mm
    chamber_length_cm = 15.1  # Length of the ionization chamber in cm
    ionization_current_pA = 2.0  # Measured current in picoamperes (pA)
    air_density_mg_cm3 = 1.293  # Density of air in mg/cm^3
    exposure_time_s = 0.02  # Exposure time for a point on the surface in seconds

    # Physical constant: W_air/e is the average energy required to create an
    # ion pair in air, per unit charge (in Joules per Coulomb).
    W_air_over_e_J_per_C = 33.97

    # --- Unit Conversions for consistency ---
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0
    ionization_current_A = ionization_current_pA * 1e-12  # Convert pA to A (C/s)
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0 # Convert mg/cm^3 to g/cm^3

    # --- Step 2: Calculate the mass of the irradiated air ---
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm
    mass_air_g = irradiated_volume_cm3 * air_density_g_cm3
    mass_air_kg = mass_air_g / 1000.0 # Convert mass to kg for dose calculation

    # --- Step 3: Calculate the dose rate in Gray/second (Gy/s) ---
    # Dose Rate = (Energy deposited per second) / mass
    # Energy deposited per second = Current (C/s) * Energy per charge (J/C)
    energy_rate_J_per_s = ionization_current_A * W_air_over_e_J_per_C
    dose_rate_Gy_per_s = energy_rate_J_per_s / mass_air_kg

    # --- Step 4: Calculate the cumulative dose in Gray (Gy) ---
    # Cumulative dose = Dose Rate * Exposure Time
    cumulative_dose_Gy = dose_rate_Gy_per_s * exposure_time_s

    # --- Step 5: Convert final results to more readable units (mGy) and print ---
    cumulative_dose_mGy = cumulative_dose_Gy * 1000.0
    dose_rate_mGy_per_s = dose_rate_Gy_per_s * 1000.0

    print("--- Calculation Breakdown ---")
    print(f"1. Mass of Irradiated Air (m):")
    print(f"   m = (Beam Area * Chamber Length) * Air Density")
    print(f"   m = (({beam_width_h_cm:.2f} cm * {beam_width_v_cm:.1f} cm) * {chamber_length_cm:.1f} cm) * {air_density_g_cm3 * 1000:.3f} mg/cm^3")
    print(f"   m = {mass_air_kg:.3e} kg")
    print("-" * 30)
    print(f"2. Dose Rate (D_rate):")
    print(f"   D_rate = (Ionization Current / Mass) * (W_air/e)")
    print(f"   D_rate = ({ionization_current_A:.1e} C/s / {mass_air_kg:.3e} kg) * {W_air_over_e_J_per_C:.2f} J/C")
    print(f"   D_rate = {dose_rate_mGy_per_s:.4f} mGy/s")
    print("-" * 30)
    print(f"3. Cumulative Surface Dose (D_cumulative):")
    print(f"   D_cumulative = Dose Rate * Exposure Time")
    print(f"   D_cumulative = {dose_rate_mGy_per_s:.4f} mGy/s * {exposure_time_s:.2f} s")
    print(f"   Final Cumulative Dose = {cumulative_dose_mGy:.6f} mGy")

if __name__ == '__main__':
    calculate_surface_dose()