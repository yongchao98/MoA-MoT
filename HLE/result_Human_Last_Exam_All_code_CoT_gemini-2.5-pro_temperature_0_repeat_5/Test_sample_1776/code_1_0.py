import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- 1. Define constants and initial values ---
    A_injection = 8.0  # mCi, desired activity at injection time
    A_cal = 10.0       # mCi, calibration activity of the vial
    V_cal = 10.0       # mL, volume of the vial
    T_half_days = 2.8  # days, half-life of Indium-111

    # Timestamps from the problem
    # The year is specified as 2024 in the problem description.
    time_injection = datetime(2024, 12, 23, 8, 0)
    time_compounding = datetime(2024, 12, 23, 4, 4)
    time_calibration = datetime(2024, 12, 26, 12, 0)

    # --- 2. Calculate the decay constant (lambda) ---
    # The formula for the decay constant is lambda = ln(2) / T_half
    lambda_per_day = math.log(2) / T_half_days

    # --- 3. Calculate activity needed at compounding time ---
    # This accounts for the decay that will occur between compounding and injection.
    # First, find the time difference t1 (from compounding to injection).
    t1_delta = time_injection - time_compounding
    t1_days = t1_delta.total_seconds() / (24 * 3600)

    # Calculate the activity needed at compounding time (A_needed).
    # We need to find the starting activity, so we use A_needed = A_final * e^(lambda * t)
    A_needed = A_injection * math.exp(lambda_per_day * t1_days)

    # --- 4. Calculate the vial's concentration at compounding time ---
    # First, find the time difference t2 (from compounding to calibration).
    t2_delta = time_calibration - time_compounding
    t2_days = t2_delta.total_seconds() / (24 * 3600)

    # Calculate the activity in the vial at compounding time (A_vial_now).
    # Compounding occurs before calibration, so the activity is higher. We use e^(lambda * t).
    A_vial_now = A_cal * math.exp(lambda_per_day * t2_days)

    # Calculate the concentration at compounding time (C_now).
    C_now = A_vial_now / V_cal

    # --- 5. Calculate the final volume to draw ---
    # Volume = Activity Needed / Concentration
    V_needed = A_needed / C_now

    # --- 6. Print the results and the final equation ---
    print("Step 1: Calculate the activity required at the time of compounding (4:04 am).")
    print(f"The dose of {A_injection:.1f} mCi is needed at {time_injection.strftime('%I:%M %p')}.")
    print(f"Time for decay between compounding and injection: {t1_delta} ({t1_days:.4f} days).")
    print(f"Required Activity = {A_injection:.1f} mCi * e^(ln(2)/{T_half_days} * {t1_days:.4f} days)")
    print(f"Required Activity = {A_needed:.4f} mCi\n")

    print("Step 2: Calculate the concentration of the Indium-111 vial at compounding time.")
    print(f"The vial is calibrated for {A_cal:.1f} mCi in {V_cal:.1f} mL at {time_calibration.strftime('%I:%M %p on %b %d, %Y')}.")
    print(f"Time from compounding to calibration: {t2_delta} ({t2_days:.4f} days).")
    print(f"Vial Activity = {A_cal:.1f} mCi * e^(ln(2)/{T_half_days} * {t2_days:.4f} days)")
    print(f"Vial Activity at compounding = {A_vial_now:.4f} mCi")
    print(f"Vial Concentration = {A_vial_now:.4f} mCi / {V_cal:.1f} mL = {C_now:.4f} mCi/mL\n")

    print("Step 3: Calculate the final volume to draw from the vial.")
    print("Volume to Draw = (Required Activity) / (Vial Concentration)")
    print(f"Volume to Draw = {A_needed:.4f} mCi / {C_now:.4f} mCi/mL")
    print(f"Volume to Draw = {V_needed:.3f} mL")
    
    # The final answer in the required format
    print(f"\n<<<{V_needed:.3f}>>>")

solve_radiopharmacy_calculation()