import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Solves the radiopharmacy calculation problem based on the user's request.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in hours
    T_half_In111_hours = 2.805 * 24

    # Calibration info for the vial
    A_calib = 10.0  # mCi
    V_calib = 10.0  # mL
    # Note: Assuming CST for all times as it's specified for the compounding time.
    time_calib = datetime(2024, 12, 26, 12, 0) # 12 PM

    # Compounding info
    time_compound = datetime(2024, 12, 23, 4, 4) # 4:04 AM CST

    # Dose info
    A_inject_target = 8.0  # mCi
    time_inject = datetime(2024, 12, 23, 8, 0) # 8 AM

    # --- Calculations ---

    # Calculate decay constant (lambda)
    lambda_val = math.log(2) / T_half_In111_hours

    print("This script calculates the required volume for a dose of Indium-111.\n")
    print(f"Constants used: Half-life of In-111 = {T_half_In111_hours:.2f} hours, Decay Constant (λ) = {lambda_val:.6f} h⁻¹\n")

    # Step 1: Calculate the vial's activity at compounding time
    print("--- Step 1: Calculate vial activity at compounding time (4:04 AM Dec 23) ---")
    time_diff_calib_to_compound = time_calib - time_compound
    t1_hours = time_diff_calib_to_compound.total_seconds() / 3600
    # Use A = A0 * e^(lambda*t) for pre-decay calculation
    A_vial_at_compound = A_calib * math.exp(lambda_val * t1_hours)
    
    print(f"Time from compounding to calibration: {time_diff_calib_to_compound} ({t1_hours:.2f} hours).")
    print("Equation: Activity_vial = Activity_calib * e^(λ * t)")
    print(f"Activity_vial = {A_calib:.1f} mCi * e^({lambda_val:.6f} * {t1_hours:.2f})")
    print(f"Calculated vial activity at compounding time: {A_vial_at_compound:.2f} mCi\n")

    # Calculate the vial's concentration at compounding time
    C_vial_at_compound = A_vial_at_compound / V_calib
    print("--- Step 2: Calculate vial concentration at compounding time ---")
    print("Equation: Concentration_vial = Activity_vial / Volume_vial")
    print(f"Concentration_vial = {A_vial_at_compound:.2f} mCi / {V_calib:.1f} mL")
    print(f"Calculated vial concentration at compounding time: {C_vial_at_compound:.3f} mCi/mL\n")

    # Step 2: Calculate the required dose activity at compounding time
    print("--- Step 3: Calculate activity needed at compounding to decay to 8 mCi by 8:00 AM ---")
    time_diff_inject_to_compound = time_inject - time_compound
    t2_hours = time_diff_inject_to_compound.total_seconds() / 3600
    # Use A_needed = A_target * e^(lambda*t) to find the activity needed earlier
    A_dose_at_compound = A_inject_target * math.exp(lambda_val * t2_hours)

    print(f"Time from compounding to injection: {time_diff_inject_to_compound} ({t2_hours:.2f} hours).")
    print("Equation: Activity_needed = Activity_target * e^(λ * t)")
    print(f"Activity_needed = {A_inject_target:.1f} mCi * e^({lambda_val:.6f} * {t2_hours:.2f})")
    print(f"Calculated activity needed at compounding time: {A_dose_at_compound:.2f} mCi\n")
    
    # Step 3: Calculate the volume to draw
    print("--- Step 4: Calculate the final volume to draw from the vial ---")
    V_draw = A_dose_at_compound / C_vial_at_compound
    print("Equation: Volume_to_draw = Activity_needed / Concentration_vial")
    print(f"Volume_to_draw = {A_dose_at_compound:.2f} mCi / {C_vial_at_compound:.3f} mCi/mL")
    print(f"\nThe final volume that must be drawn into the syringe is {V_draw:.2f} mL.")


solve_radiopharmacy_calculation()

# The final answer is calculated by the script. Re-running the calculation for the final tag:
T_half_In111_hours = 2.805 * 24
lambda_val = math.log(2) / T_half_In111_hours
time_calib = datetime(2024, 12, 26, 12, 0)
time_compound = datetime(2024, 12, 23, 4, 4)
time_inject = datetime(2024, 12, 23, 8, 0)
t1_hours = (time_calib - time_compound).total_seconds() / 3600
A_vial_at_compound = 10.0 * math.exp(lambda_val * t1_hours)
C_vial_at_compound = A_vial_at_compound / 10.0
t2_hours = (time_inject - time_compound).total_seconds() / 3600
A_dose_at_compound = 8.0 * math.exp(lambda_val * t2_hours)
V_draw = A_dose_at_compound / C_vial_at_compound
print(f"<<<{V_draw:.2f}>>>")