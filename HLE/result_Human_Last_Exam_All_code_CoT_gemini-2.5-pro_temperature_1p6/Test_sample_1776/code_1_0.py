import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride to be drawn for a patient dose.
    """
    # --- Input values from the problem statement ---
    calib_activity = 10.0  # mCi
    calib_volume = 10.0  # mL
    calib_datetime = datetime(2024, 12, 26, 12, 0)
    
    ordered_activity_at_injection = 8.0  # mCi
    injection_datetime = datetime(2024, 12, 23, 8, 0)
    
    compounding_datetime = datetime(2024, 12, 23, 4, 4)
    
    # Half-life of Indium-111
    half_life_in_days = 2.80
    half_life_in_hours = half_life_in_days * 24

    # --- Step 1: Calculate Time Differences in Hours ---
    print("### Step 1: Calculate Time Differences ###")
    
    # Time from compounding to calibration (for vial activity calculation)
    time_diff_pre_calib = calib_datetime - compounding_datetime
    t_pre_calib_hours = time_diff_pre_calib.total_seconds() / 3600
    print(f"Time from compounding to vial calibration: {t_pre_calib_hours:.2f} hours")

    # Time from compounding to injection (for dose decay calculation)
    time_diff_decay = injection_datetime - compounding_datetime
    t_decay_hours = time_diff_decay.total_seconds() / 3600
    print(f"Time from compounding to patient injection: {t_decay_hours:.2f} hours\n")

    # --- Step 2: Calculate Decay Constant (λ) ---
    print("### Step 2: Calculate Decay Constant (λ) for Indium-111 ###")
    decay_constant = math.log(2) / half_life_in_hours
    print(f"The half-life of In-111 is {half_life_in_hours} hours.")
    print(f"Decay constant λ = ln(2) / T_half = {math.log(2):.4f} / {half_life_in_hours} h = {decay_constant:.6f} per hour\n")
    
    # --- Step 3: Calculate the Required Activity at Compounding Time ---
    # This is the activity needed so it decays to the ordered dose at injection time.
    # Formula: A_needed = A_ordered * e^(λ * t_decay)
    print("### Step 3: Calculate Activity Needed at Compounding Time ###")
    activity_needed_at_compounding = ordered_activity_at_injection * math.exp(decay_constant * t_decay_hours)
    print("Activity Needed = Ordered Activity * e^(λ * t_decay)")
    print(f"Activity Needed = {ordered_activity_at_injection} mCi * e^({decay_constant:.6f} * {t_decay_hours:.2f})")
    print(f"Result: {activity_needed_at_compounding:.4f} mCi\n")
    
    # --- Step 4: Calculate the Vial's Actual Activity at Compounding Time ---
    # The vial is calibrated for a future date, so its activity now is higher.
    # Formula: A_vial_now = A_calibrated * e^(λ * t_pre_calib)
    print("### Step 4: Calculate Vial's Activity at Compounding Time ###")
    activity_in_vial_at_compounding = calib_activity * math.exp(decay_constant * t_pre_calib_hours)
    print("Vial Activity = Calibrated Activity * e^(λ * t_pre_calib)")
    print(f"Vial Activity = {calib_activity} mCi * e^({decay_constant:.6f} * {t_pre_calib_hours:.2f})")
    print(f"Result: {activity_in_vial_at_compounding:.4f} mCi\n")

    # --- Step 5: Calculate the Vial's Concentration at Compounding Time ---
    print("### Step 5: Calculate Vial's Concentration at Compounding Time ###")
    concentration_at_compounding = activity_in_vial_at_compounding / calib_volume
    print("Vial Concentration = Vial Activity / Vial Volume")
    print(f"Vial Concentration = {activity_in_vial_at_compounding:.4f} mCi / {calib_volume} mL")
    print(f"Result: {concentration_at_compounding:.4f} mCi/mL\n")
    
    # --- Step 6: Calculate the Final Volume to Draw ---
    print("### Step 6: Calculate Final Volume to Draw ###")
    volume_to_draw = activity_needed_at_compounding / concentration_at_compounding
    print("Volume to Draw = Activity Needed at Compounding / Vial Concentration at Compounding")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL\n")
    
    print("---------------------------------------------------------------")
    print(f"Final Answer:")
    print(f"The volume of Indium 111 chloride that must be drawn is {volume_to_draw:.2f} mL.")
    print("---------------------------------------------------------------")
    
    return volume_to_draw

# Execute the function
if __name__ == '__main__':
    final_volume = solve_radiopharmacy_calculation()
    # The final numerical result for the platform
    # print(f"<<<{final_volume:.2f}>>>")

<<<3.65>>>