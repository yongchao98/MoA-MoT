import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride required for a patient dose.
    """
    # --- Step 1: Define initial parameters and constants ---
    
    # Dose Information
    A_ordered = 8.0  # mCi
    t_inject_str = "2024-12-23 08:00:00"
    t_compound_str = "2024-12-23 04:04:00"

    # Stock Vial Information
    A_calib = 10.0  # mCi
    V_calib = 10.0  # mL
    t_calib_str = "2024-12-26 12:00:00"

    # Radionuclide Information (Indium-111)
    T_half_days = 2.80  # days

    # Convert times to datetime objects for accurate calculations
    t_inject = datetime.datetime.strptime(t_inject_str, "%Y-%m-%d %H:%M:%S")
    t_compound = datetime.datetime.strptime(t_compound_str, "%Y-%m-%d %H:%M:%S")
    t_calib = datetime.datetime.strptime(t_calib_str, "%Y-%m-%d %H:%M:%S")
    
    # Calculate half-life in hours and the decay constant (lambda)
    T_half_h = T_half_days * 24.0
    lambda_h = math.log(2) / T_half_h
    
    print("--- Plan ---")
    print("1. Calculate the activity needed at the time of compounding (4:04 am).")
    print("2. Calculate the activity and concentration of the stock vial at the same time.")
    print("3. Divide the needed activity by the vial's concentration to find the required volume.")
    print("\n--- Calculations ---")
    
    # --- Step 2: Calculate activity needed at compounding time ---
    
    # Time from compounding to injection
    dt_kit = t_inject - t_compound
    dt_kit_h = dt_kit.total_seconds() / 3600.0
    
    # Activity needed at compounding must be higher to account for decay before injection
    # Formula: A_initial = A_final * e^(lambda * t)
    A_needed_at_compounding = A_ordered * math.exp(lambda_h * dt_kit_h)
    
    print(f"1. Dose required at injection time (8:00 am): {A_ordered:.2f} mCi")
    print(f"   Time from compounding to injection: {dt_kit_h:.2f} hours")
    print(f"   Activity needed at compounding time (4:04 am): {A_needed_at_compounding:.2f} mCi")
    print("-" * 20)

    # --- Step 3: Calculate the vial's concentration at compounding time ---

    # Time from compounding to vial calibration
    dt_vial = t_calib - t_compound
    dt_vial_h = dt_vial.total_seconds() / 3600.0
    
    # Activity of the vial at compounding time (it's higher because it's before calibration)
    # Formula: A_initial = A_final * e^(lambda * t)
    A_vial_at_compounding = A_calib * math.exp(lambda_h * dt_vial_h)
    
    # Concentration of the vial at compounding time
    C_vial_at_compounding = A_vial_at_compounding / V_calib
    
    print(f"2. Vial calibrated activity at {t_calib.strftime('%b %d, %I:%M %p')}: {A_calib:.2f} mCi in {V_calib:.1f} mL")
    print(f"   Time from compounding to calibration: {dt_vial_h:.2f} hours")
    print(f"   Vial activity at compounding time: {A_vial_at_compounding:.2f} mCi")
    print(f"   Vial concentration at compounding time: {C_vial_at_compounding:.2f} mCi/mL")
    print("-" * 20)
    
    # --- Step 4: Calculate the final volume to draw ---
    
    volume_to_draw = A_needed_at_compounding / C_vial_at_compounding

    print("3. To find the volume to draw, divide the needed activity by the vial's concentration:")
    print("\n--- Final Equation ---")
    print(f"Volume to Draw (mL) = Activity Needed (mCi) / Vial Concentration (mCi/mL)")
    print(f"Volume = {A_needed_at_compounding:.2f} mCi / {C_vial_at_compounding:.2f} mCi/mL")
    print(f"\nFinal Volume to Draw: {volume_to_draw:.2f} mL")
    
    return volume_to_draw

# Run the calculation and print the final result in the requested format
final_volume = solve_radiopharmacy_calculation()
# The final result is returned for format compliance, but already printed above.
# print(f"<<<{final_volume:.2f}>>>")