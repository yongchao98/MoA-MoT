import math
import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 for a patient dose,
    accounting for radioactive decay.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in hours
    T_half_hours = 67.32

    # Desired dose for the patient in mCi
    desired_dose_mci = 8.0

    # Calibration information
    calibrated_activity_mci = 10.0
    calibrated_volume_ml = 10.0
    calibration_datetime = datetime.datetime(2024, 12, 26, 12, 0) # 12 pm on Dec 26

    # Compounding time
    compounding_datetime = datetime.datetime(2024, 12, 23, 4, 4) # 4:04 am on Dec 23

    # --- Step 1: Calculate the time difference in hours ---
    time_difference = calibration_datetime - compounding_datetime
    time_diff_hours = time_difference.total_seconds() / 3600.0

    # --- Step 2: Calculate the decay constant (lambda) ---
    # lambda = ln(2) / T_half
    decay_constant = math.log(2) / T_half_hours

    # --- Step 3: Calculate the activity at compounding time ---
    # Formula for pre-calibration: A_compound = A_cal * e^(lambda * t)
    activity_at_compounding = calibrated_activity_mci * math.exp(decay_constant * time_diff_hours)

    # --- Step 4: Calculate the concentration at compounding time ---
    # Concentration = Total Activity / Total Volume
    concentration_at_compounding = activity_at_compounding / calibrated_volume_ml

    # --- Step 5: Calculate the final volume needed ---
    # Volume = Desired Dose / Concentration
    volume_needed_ml = desired_dose_mci / concentration_at_compounding

    # --- Final Output ---
    print("To find the required volume, we use the formula for radioactive decay.")
    print("First, we calculate the activity concentration (mCi/mL) at the time of compounding.")
    print("\nEquation for concentration at compounding time:")
    print("Concentration = (Calibrated Activity / Calibrated Volume) * e^((ln(2) / Half-Life) * Time Difference)")
    print(f"Concentration = ({calibrated_activity_mci} mCi / {calibrated_volume_ml} mL) * e^((ln(2) / {T_half_hours} hr) * {time_diff_hours:.3f} hr)")
    
    decay_factor = math.exp(decay_constant * time_diff_hours)
    print(f"Concentration = {calibrated_activity_mci / calibrated_volume_ml} mCi/mL * e^({decay_constant * time_diff_hours:.4f})")
    print(f"Concentration = {calibrated_activity_mci / calibrated_volume_ml} mCi/mL * {decay_factor:.4f}")
    print(f"Concentration = {concentration_at_compounding:.4f} mCi/mL")

    print("\nNext, we calculate the volume needed to get the desired dose:")
    print("Volume to Draw = Desired Dose / Concentration at Compounding Time")
    print(f"Volume to Draw = {desired_dose_mci} mCi / {concentration_at_compounding:.4f} mCi/mL")
    print(f"Volume to Draw = {volume_needed_ml:.3f} mL")


solve_radiopharmacy_calculation()
<<<3.512>>>