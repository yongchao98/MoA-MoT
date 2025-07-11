import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a dose.
    """
    # --- Given Information ---
    # Radionuclide properties
    half_life_days = 2.80
    half_life_hours = half_life_days * 24

    # Dose information
    dose_activity_mci = 8.0
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Kit (stock vial) information
    cal_activity_mci = 10.0
    cal_volume_ml = 10.0
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Compounding time
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # --- Calculations ---

    # Decay constant (lambda)
    decay_constant_per_hour = math.log(2) / half_life_hours

    print("--- Step 1: Calculate the activity needed at compounding time ---")
    # Time difference for dose decay (from compounding to injection)
    dt_dose = injection_time - compounding_time
    dt_dose_hours = dt_dose.total_seconds() / 3600
    
    # Pre-cay calculation: find activity needed at the earlier compounding time
    activity_needed_at_compounding = dose_activity_mci * math.exp(decay_constant_per_hour * dt_dose_hours)
    
    print("The dose of 8.00 mCi is needed at 8:00 am.")
    print("To account for decay between compounding (4:04 am) and injection (8:00 am), we must draw more activity.")
    print(f"Time for decay (t1): {dt_dose_hours:.3f} hours")
    print("Equation: Activity_Needed = Dose_Activity * e^(lambda * t1)")
    print(f"Activity_Needed = {dose_activity_mci:.2f} mCi * e^({decay_constant_per_hour:.6f} * {dt_dose_hours:.3f}) = {activity_needed_at_compounding:.3f} mCi")
    print("-" * 20)

    print("\n--- Step 2: Calculate the stock vial concentration at compounding time ---")
    # Time difference for vial decay (from compounding to calibration)
    dt_vial = calibration_time - compounding_time
    dt_vial_hours = dt_vial.total_seconds() / 3600

    # Pre-cay calculation: find the vial's activity at the earlier compounding time
    vial_activity_at_compounding = cal_activity_mci * math.exp(decay_constant_per_hour * dt_vial_hours)
    vial_concentration_at_compounding = vial_activity_at_compounding / cal_volume_ml

    print("The stock vial is calibrated for 10.0 mCi at a future date (Dec 26, 12:00 pm).")
    print("Its activity at the earlier compounding time (Dec 23, 4:04 am) will be higher.")
    print(f"Time from compounding to calibration (t2): {dt_vial_hours:.3f} hours")
    print("Equation: Vial_Activity = Cal_Activity * e^(lambda * t2)")
    print(f"Vial_Activity = {cal_activity_mci:.2f} mCi * e^({decay_constant_per_hour:.6f} * {dt_vial_hours:.3f}) = {vial_activity_at_compounding:.3f} mCi")
    print("\nEquation: Vial_Concentration = Vial_Activity / Vial_Volume")
    print(f"Vial_Concentration = {vial_activity_at_compounding:.3f} mCi / {cal_volume_ml:.1f} mL = {vial_concentration_at_compounding:.4f} mCi/mL")
    print("-" * 20)

    print("\n--- Step 3: Calculate the final volume to draw ---")
    volume_to_draw_ml = activity_needed_at_compounding / vial_concentration_at_compounding

    print("Equation: Volume_to_Draw = Activity_Needed / Vial_Concentration")
    print(f"Volume_to_Draw = {activity_needed_at_compounding:.3f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")
    print("-" * 20)

    print(f"\nFinal Answer: The volume of Indium 111 chloride that must be drawn is {volume_to_draw_ml:.2f} mL.")

solve_radiopharmacy_calculation()
<<<3.65>>>