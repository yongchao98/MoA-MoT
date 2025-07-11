import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Input Parameters ---
    
    # Half-life of Indium-111 is 2.805 days
    T_HALF_DAYS = 2.805
    T_HALF_HOURS = T_HALF_DAYS * 24

    # Vial calibration details
    CALIBRATION_ACTIVITY_mCi = 10.0
    CALIBRATION_VOLUME_mL = 10.0
    time_calibration = datetime.datetime(2024, 12, 26, 12, 0) # 12 pm on Dec 26

    # Compounding and injection details
    time_compounding = datetime.datetime(2024, 12, 23, 4, 4)  # 4:04 am on Dec 23
    time_injection = datetime.datetime(2024, 12, 23, 8, 0)   # 8 am on Dec 23
    DOSE_AT_INJECTION_mCi = 8.0

    # --- Calculations ---
    
    # The decay constant (lambda)
    DECAY_CONSTANT_per_hour = math.log(2) / T_HALF_HOURS

    # Step 1: Find the vial's concentration at the time of compounding.
    # Because compounding is before calibration, the activity will be higher.
    # The formula is A = A0 * e^(lambda * t), where A0 is the future calibrated activity.
    time_from_compounding_to_cal = time_calibration - time_compounding
    hours_from_compounding_to_cal = time_from_compounding_to_cal.total_seconds() / 3600
    
    activity_vial_at_compounding = CALIBRATION_ACTIVITY_mCi * math.exp(DECAY_CONSTANT_per_hour * hours_from_compounding_to_cal)
    concentration_at_compounding = activity_vial_at_compounding / CALIBRATION_VOLUME_mL

    # Step 2: Find the required activity at compounding time for the patient dose.
    # We need to draw a higher activity to account for decay between compounding and injection.
    time_from_compounding_to_injection = time_injection - time_compounding
    hours_from_compounding_to_injection = time_from_compounding_to_injection.total_seconds() / 3600

    activity_needed_at_compounding = DOSE_AT_INJECTION_mCi * math.exp(DECAY_CONSTANT_per_hour * hours_from_compounding_to_injection)

    # Step 3: Calculate the volume to draw.
    # Volume = Activity Needed / Concentration Available
    volume_to_draw_mL = activity_needed_at_compounding / concentration_at_compounding

    # --- Output Results ---
    print("--- Radiopharmacy Calculation ---")
    print(f"Indium-111 Half-life: {T_HALF_HOURS:.2f} hours")
    print("-" * 35)

    print("Step 1: Calculate vial concentration at compounding time (4:04 am Dec 23)")
    print(f"Time from compounding to calibration: {hours_from_compounding_to_cal:.2f} hours")
    print(f"Vial activity at compounding: {activity_vial_at_compounding:.4f} mCi in {CALIBRATION_VOLUME_mL} mL")
    print(f"Vial concentration at compounding: {concentration_at_compounding:.4f} mCi/mL")
    print("-" * 35)

    print("Step 2: Calculate required dose activity at compounding time")
    print(f"Time from compounding to injection: {hours_from_compounding_to_injection:.2f} hours")
    print(f"Required activity to have {DOSE_AT_INJECTION_mCi} mCi at 8 am: {activity_needed_at_compounding:.4f} mCi")
    print("-" * 35)

    print("Step 3: Calculate the final volume to draw")
    print("Formula: Volume = (Required Activity) / (Vial Concentration)")
    print(f"Volume = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"\nFinal volume to draw: {volume_to_draw_mL:.2f} mL")


solve_radiopharmacy_calculation()
<<<3.65>>>