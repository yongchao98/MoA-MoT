import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for an Octreascan kit.
    The calculation proceeds in three main steps:
    1. Determine the activity needed at compounding time to ensure the dose is correct at injection time.
    2. Determine the concentration of the stock vial at compounding time based on its future calibration point.
    3. Divide the needed activity by the stock concentration to find the required volume.
    """
    # --- Input Parameters ---
    T_half_days = 2.805  # Half-life of Indium-111 in days
    desired_activity_at_injection = 8.0  # mCi
    injection_time_str = "2024-12-23 08:00:00"
    initial_activity_at_calibration = 10.0  # mCi
    initial_volume = 10.0  # mL
    calibration_time_str = "2024-12-26 12:00:00"
    compounding_time_str = "2024-12-23 04:04:00"

    # --- Calculations ---

    # Convert T_half to hours and calculate the decay constant (lambda)
    T_half_hours = T_half_days * 24
    decay_constant = math.log(2) / T_half_hours  # per hour

    # Parse all time strings into datetime objects
    injection_time = datetime.strptime(injection_time_str, "%Y-%m-%d %H:%M:%S")
    calibration_time = datetime.strptime(calibration_time_str, "%Y-%m-%d %H:%M:%S")
    compounding_time = datetime.strptime(compounding_time_str, "%Y-%m-%d %H:%M:%S")

    # Time differences in hours
    time_decay_to_injection = (injection_time - compounding_time).total_seconds() / 3600
    time_from_compounding_to_calibration = (calibration_time - compounding_time).total_seconds() / 3600

    # Step 1: Calculate activity needed at compounding time
    # This accounts for decay between compounding and injection. Formula: A_needed = A_final * exp(lambda * t)
    activity_needed_at_compounding = desired_activity_at_injection * math.exp(decay_constant * time_decay_to_injection)

    # Step 2: Calculate stock vial activity and concentration at compounding time
    # This accounts for the vial being calibrated at a future date. Formula: A_past = A_future * exp(lambda * t)
    vial_activity_at_compounding = initial_activity_at_calibration * math.exp(decay_constant * time_from_compounding_to_calibration)
    vial_concentration_at_compounding = vial_activity_at_compounding / initial_volume  # mCi/mL

    # Step 3: Calculate the final volume to draw
    volume_to_draw = activity_needed_at_compounding / vial_concentration_at_compounding

    # --- Output ---
    print("This calculation determines the volume of Indium-111 to draw for the patient's dose.")

    print("\n--- Step 1: Calculate Required Activity at Compounding Time (04:04 am) ---")
    print(f"The dose of {desired_activity_at_injection:.1f} mCi is needed at {injection_time.strftime('%I:%M %p')}.")
    print(f"To account for decay over {time_decay_to_injection:.2f} hours until injection, a higher activity is needed at compounding time.")
    print(f"Activity Needed = Desired Activity * e^(λ * t)")
    print(f"Activity Needed = {desired_activity_at_injection:.1f} mCi * e^({decay_constant:.6f} * {time_decay_to_injection:.2f} h)")
    print(f"Activity Needed = {activity_needed_at_compounding:.4f} mCi")

    print("\n--- Step 2: Calculate Stock Vial Concentration at Compounding Time (04:04 am) ---")
    print(f"The stock vial is calibrated for {initial_activity_at_calibration:.1f} mCi at {calibration_time.strftime('%I:%M %p on %b %d')}.")
    print(f"The time from compounding to calibration is {time_from_compounding_to_calibration:.2f} hours.")
    print("The vial's activity at the earlier compounding time must be calculated.")
    print(f"Vial Activity = Calibration Activity * e^(λ * t)")
    print(f"Vial Activity = {initial_activity_at_calibration:.1f} mCi * e^({decay_constant:.6f} * {time_from_compounding_to_calibration:.2f} h)")
    print(f"Vial Activity at Compounding = {vial_activity_at_compounding:.4f} mCi")
    print(f"Vial Concentration = {vial_activity_at_compounding:.4f} mCi / {initial_volume:.1f} mL = {vial_concentration_at_compounding:.4f} mCi/mL")

    print("\n--- Step 3: Calculate Final Volume to Draw ---")
    print("The volume to draw is the required activity divided by the vial's concentration.")
    print(f"Volume = Activity Needed / Vial Concentration")
    print(f"Volume = {activity_needed_at_compounding:.4f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")
    print(f"\nFinal Volume to be drawn: {volume_to_draw:.2f} mL")


solve_radiopharmacy_calculation()