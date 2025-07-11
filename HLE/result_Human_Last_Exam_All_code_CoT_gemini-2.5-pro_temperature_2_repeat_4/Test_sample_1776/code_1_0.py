import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a dose,
    accounting for radioactive decay.
    """
    # Step 1: Define constants and initial values
    T_HALF_HOURS = 2.805 * 24  # Half-life of Indium-111 in hours (67.32 hours)
    DECAY_CONSTANT = math.log(2) / T_HALF_HOURS  # Decay constant (lambda)

    # Dose information
    INJECTION_DOSE_ACTIVITY = 8.0  # mCi
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Vial information
    CALIBRATION_ACTIVITY = 10.0  # mCi
    CALIBRATION_VOLUME = 10.0  # mL
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Compounding time
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # --- Print initial values ---
    print("--- Radiopharmacy Calculation for Indium-111 ---")
    print("\nStep 1: Constants and Initial Values")
    print(f"Half-life of Indium-111 (T½): {T_HALF_HOURS:.2f} hours")
    print(f"Decay constant (λ): {DECAY_CONSTANT:.6f} h⁻¹")
    print(f"Dose required at injection: {INJECTION_DOSE_ACTIVITY} mCi at {injection_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Vial calibration: {CALIBRATION_ACTIVITY} mCi in {CALIBRATION_VOLUME} mL at {calibration_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Compounding time: {compounding_time.strftime('%Y-%m-%d %H:%M')}")

    # Step 2: Calculate Time Differences
    time_from_compounding_to_calibration = calibration_time - compounding_time
    t1_hours = time_from_compounding_to_calibration.total_seconds() / 3600

    time_from_compounding_to_injection = injection_time - compounding_time
    t2_hours = time_from_compounding_to_injection.total_seconds() / 3600
    
    print("\nStep 2: Time Calculations")
    print(f"Time from compounding to vial calibration (t1): {t1_hours:.2f} hours")
    print(f"Time from compounding to dose injection (t2): {t2_hours:.2f} hours")

    # Step 3: Calculate Vial Activity at Compounding Time
    # Formula: A_at_compounding = A_at_calibration * e^(lambda * t1)
    # This accounts for the decay that will happen between compounding and calibration.
    vial_activity_at_compounding = CALIBRATION_ACTIVITY * math.exp(DECAY_CONSTANT * t1_hours)
    
    print("\nStep 3: Calculate Vial Activity at Compounding Time")
    print(f"This uses the formula: A_compounding = A_calibration * e^(λ * t1)")
    print(f"Vial Activity = {CALIBRATION_ACTIVITY} mCi * e^({DECAY_CONSTANT:.6f} * {t1_hours:.2f})")
    print(f"Vial Activity at Compounding Time: {vial_activity_at_compounding:.2f} mCi")

    # Step 4: Calculate Activity Needed for Dose at Compounding Time
    # Formula: A_to_draw = A_at_injection * e^(lambda * t2)
    # This ensures that after decay (t2), the dose has the correct activity.
    activity_to_draw = INJECTION_DOSE_ACTIVITY * math.exp(DECAY_CONSTANT * t2_hours)
    
    print("\nStep 4: Calculate Activity Needed for Dose at Compounding Time")
    print(f"This uses the formula: A_draw = A_injection * e^(λ * t2)")
    print(f"Activity to Draw = {INJECTION_DOSE_ACTIVITY} mCi * e^({DECAY_CONSTANT:.6f} * {t2_hours:.2f})")
    print(f"Activity Needed at Compounding Time: {activity_to_draw:.2f} mCi")

    # Step 5: Calculate Volume to Draw
    concentration_at_compounding = vial_activity_at_compounding / CALIBRATION_VOLUME
    volume_to_draw = activity_to_draw / concentration_at_compounding

    print("\nStep 5: Calculate Final Volume to Draw")
    print(f"Vial Concentration at Compounding = {vial_activity_at_compounding:.2f} mCi / {CALIBRATION_VOLUME} mL = {concentration_at_compounding:.3f} mCi/mL")
    print("\n--- Final Answer ---")
    
    # Print the final equation with all the calculated numbers
    print("Final Equation:")
    print(f"Volume to Draw (mL) = [Activity Needed at Compounding (mCi)] / [Vial Concentration at Compounding (mCi/mL)]")
    print(f"Volume to Draw (mL) = {activity_to_draw:.2f} / ({vial_activity_at_compounding:.2f} / {CALIBRATION_VOLUME})")
    print(f"The required volume to draw is: {volume_to_draw:.2f} mL")
    
    return volume_to_draw

if __name__ == '__main__':
    final_volume = solve_radiopharmacy_calculation()
    print(f"\n<<<{final_volume:.2f}>>>")