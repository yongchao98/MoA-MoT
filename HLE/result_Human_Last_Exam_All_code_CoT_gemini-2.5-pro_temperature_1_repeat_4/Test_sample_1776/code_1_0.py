import math
import datetime

def solve_nuclear_pharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride to be drawn for a patient dose.
    """
    # --- Given Information ---
    # Dose Information
    ordered_activity = 8.0  # mCi
    injection_time = datetime.datetime(2024, 12, 23, 8, 0)

    # Vial Information
    calibration_activity = 10.0  # mCi
    calibration_volume = 10.0  # mL
    calibration_time = datetime.datetime(2024, 12, 26, 12, 0)

    # Compounding Information
    compounding_time = datetime.datetime(2024, 12, 23, 4, 4)

    # Radionuclide Constant
    half_life_days = 2.8
    half_life_hours = half_life_days * 24

    # --- Calculations ---
    print("Step 1: Calculate the activity needed for the dose at compounding time.")
    # Decay constant (lambda)
    decay_constant_hr = 0.693 / half_life_hours
    
    # Time difference between compounding and injection
    time_diff_dose = injection_time - compounding_time
    time_diff_dose_hours = time_diff_dose.total_seconds() / 3600
    
    # Activity needed at compounding (pre-calibrating the dose)
    # A_0 = A * e^(lambda * t)
    activity_needed_at_compounding = ordered_activity * math.exp(decay_constant_hr * time_diff_dose_hours)
    
    print(f"Time from compounding (4:04 am) to injection (8:00 am): {time_diff_dose_hours:.2f} hours")
    print(f"The dose must have {activity_needed_at_compounding:.4f} mCi at compounding time to decay to {ordered_activity} mCi at injection time.\n")

    print("Step 2: Calculate the vial's concentration at compounding time.")
    # Time difference between compounding and calibration
    time_diff_vial = calibration_time - compounding_time
    time_diff_vial_hours = time_diff_vial.total_seconds() / 3600
    
    # Activity of the entire vial at compounding time
    vial_activity_at_compounding = calibration_activity * math.exp(decay_constant_hr * time_diff_vial_hours)
    
    # Concentration of the vial at compounding time
    vial_concentration_at_compounding = vial_activity_at_compounding / calibration_volume
    
    print(f"Time from compounding (4:04 am Dec 23) to calibration (12:00 pm Dec 26): {time_diff_vial_hours:.2f} hours")
    print(f"The vial's total activity at compounding time is {vial_activity_at_compounding:.4f} mCi.")
    print(f"The vial's concentration at compounding time is {vial_concentration_at_compounding:.4f} mCi/mL.\n")

    print("Step 3: Calculate the final volume to draw from the vial.")
    # Volume = Required Activity / Concentration
    volume_to_draw = activity_needed_at_compounding / vial_concentration_at_compounding
    
    print("The final calculation is:")
    print(f"Volume = (Activity needed for dose) / (Vial Concentration)")
    print(f"Volume = {activity_needed_at_compounding:.4f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")
    print(f"The required volume to draw is {volume_to_draw:.2f} mL.")
    
    # Return the final numerical answer for the platform.
    return round(volume_to_draw, 2)

# Execute the function and print the final result in the required format
final_answer = solve_nuclear_pharmacy_calculation()
# The problem asks for the final answer in a specific format at the end.
# The text is printed by the function above.
print(f"\n<<< {final_answer} >>>")
