import sys

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ion chamber measurements.
    """

    # --- Given Parameters ---
    # Convert all initial values to SI units for consistency (meters, kilograms, seconds, Coulombs, etc.)

    # Photon beam size at focus (horizontal by vertical)
    beam_h_focal_m = 0.3 / 1000  # 0.3 mm to meters
    beam_v_focal_m = 6 / 1000    # 6 mm to meters

    # Ionization chamber properties
    chamber_length_m = 15.1 / 100 # 15.1 cm to meters
    # Density of air: 1.293 mg/cm^3 -> kg/m^3
    # 1.293 * (1e-6 kg / 1 mg) * (1e6 cm^3 / 1 m^3) = 1.293 kg/m^3
    air_density_kg_m3 = 1.293

    # Measurement readings
    ion_current_A = 2.0e-12  # 2.0 pA to Amperes (C/s)
    exposure_time_s = 0.02   # in seconds

    # --- Physical Constants ---
    # Average energy to create an ion pair in air (W_air / e)
    w_over_e_air_J_C = 33.97  # in Joules/Coulomb

    # --- Step 1: Calculate the irradiated volume of air ---
    beam_area_m2 = beam_h_focal_m * beam_v_focal_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m

    # --- Step 2: Calculate the mass of the irradiated air ---
    irradiated_mass_kg = irradiated_volume_m3 * air_density_kg_m3

    # --- Step 3 & 4: Calculate the dose rate in air (and thus tissue) ---
    # Dose Rate (Gy/s) = (Energy per second) / mass = ( (Charge per second) * (Energy per charge) ) / mass
    # Dose Rate = (Current * (W/e)) / mass
    dose_rate_Gy_s = (ion_current_A / irradiated_mass_kg) * w_over_e_air_J_C

    # --- Step 5: Calculate the cumulative dose ---
    # Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

    # --- Final Output ---
    # The user requested to show each number in the final equation.
    # Let's show the calculation for dose rate first, then the cumulative dose.
    
    # Using f-strings to format scientific notation for clarity.
    # We round intermediate results for display purposes only.
    # The actual calculation uses the full precision numbers.
    
    print("Step 1: Calculate Mass of Irradiated Air")
    print(f"Mass = (Beam Area [m^2] * Chamber Length [m]) * Air Density [kg/m^3]")
    print(f"Mass = (({beam_h_focal_m:.1e} * {beam_v_focal_m:.1e}) * {chamber_length_m:.3f}) * {air_density_kg_m3:.4f}")
    print(f"Mass = {irradiated_mass_kg:.4e} kg\n")

    print("Step 2: Calculate Dose Rate")
    print(f"Dose Rate [Gy/s] = (Ionization Current [A] / Mass [kg]) * (W/e for air [J/C])")
    print(f"Dose Rate = ({ion_current_A:.1e} A / {irradiated_mass_kg:.4e} kg) * {w_over_e_air_J_C:.2f} J/C")
    print(f"Dose Rate = {dose_rate_Gy_s:.4e} Gy/s\n")

    print("Step 3: Calculate Cumulative Surface Dose")
    print(f"Cumulative Dose [Gy] = Dose Rate [Gy/s] * Exposure Time [s]")
    print(f"Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s:.2f} s")
    print(f"Cumulative Dose = {cumulative_dose_Gy:.4e} Gy\n")
    
    # Store the final numerical answer in the required format for the calling system.
    # Do not print this special formatted string to the console for the user.
    # It will be extracted from the raw response.
    sys.stdout.flush() # ensure previous prints are out
    final_answer_string_for_extraction = f"<<<{cumulative_dose_Gy:.3e}>>>"
    # To prevent this from printing, we just define it.
    # The calling environment will see it. 
    # For a runnable script, we could write this to a file or handle differently.
    # For this interactive context, let's just add it to the final printed output.
    print("Final Answer for Submission:")
    print(f"{cumulative_dose_Gy:.3e}")


calculate_surface_dose()