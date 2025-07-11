import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of a radioactive isotope for a patient dose,
    accounting for decay from compounding time to injection time and the decay
    of the source vial relative to its calibration date.
    """

    # --- 1. Define constants and timestamps ---
    HALF_LIFE_DAYS = 2.805
    ORDERED_DOSE_mCi = 8.0
    VIAL_CALIBRATION_ACTIVITY_mCi = 10.0
    VIAL_VOLUME_mL = 10.0

    CALIBRATION_TIME = datetime(2024, 12, 26, 12, 0)   # 12:00 pm Dec 26
    INJECTION_TIME = datetime(2024, 12, 23, 8, 0)     # 08:00 am Dec 23
    COMPOUNDING_TIME = datetime(2024, 12, 23, 4, 4)    # 04:04 am Dec 23

    # Convert half-life to hours and calculate the decay constant (lambda)
    half_life_hours = HALF_LIFE_DAYS * 24
    decay_constant_lambda = math.log(2) / half_life_hours

    # --- 2. Calculate activity needed at compounding time ---
    # This is the activity required for the dose, which will decay to 8 mCi by injection time.
    # It is a "pre-decay" calculation, so we use a positive exponent.
    time_diff_dose_decay = INJECTION_TIME - COMPOUNDING_TIME
    time_diff_dose_decay_hours = time_diff_dose_decay.total_seconds() / 3600
    activity_needed_at_compounding = ORDERED_DOSE_mCi * math.exp(decay_constant_lambda * time_diff_dose_decay_hours)

    # --- 3. Calculate stock vial concentration at compounding time ---
    # The vial's activity is known at a future calibration time. Its activity now will be higher.
    time_diff_vial_decay = CALIBRATION_TIME - COMPOUNDING_TIME
    time_diff_vial_decay_hours = time_diff_vial_decay.total_seconds() / 3600
    vial_activity_at_compounding = VIAL_CALIBRATION_ACTIVITY_mCi * math.exp(decay_constant_lambda * time_diff_vial_decay_hours)
    vial_concentration_at_compounding = vial_activity_at_compounding / VIAL_VOLUME_mL

    # --- 4. Calculate the final volume to draw ---
    volume_to_draw_ml = activity_needed_at_compounding / vial_concentration_at_compounding

    # --- Print the results step-by-step ---
    print("Step-by-step Radiopharmacy Calculation:\n")
    print(f"1. Decay Constant (λ) for In-111 (Half-life {half_life_hours:.2f} hrs): {decay_constant_lambda:.6f}/hr\n")

    print("2. Activity required for the patient dose at compounding time (4:04 am):")
    print(f"   - To have {ORDERED_DOSE_mCi:.2f} mCi at 8:00 am ({time_diff_dose_decay_hours:.2f} hours later), we need:")
    print(f"   - A_needed = {ORDERED_DOSE_mCi:.2f} mCi * e^({decay_constant_lambda:.6f} * {time_diff_dose_decay_hours:.2f}h) = {activity_needed_at_compounding:.2f} mCi\n")

    print("3. Concentration of the stock vial at compounding time (4:04 am):")
    print(f"   - The vial will decay for {time_diff_vial_decay_hours:.2f} hours to reach {VIAL_CALIBRATION_ACTIVITY_mCi:.2f} mCi.")
    print(f"   - Current total activity in vial = {VIAL_CALIBRATION_ACTIVITY_mCi:.2f} mCi * e^({decay_constant_lambda:.6f} * {time_diff_vial_decay_hours:.2f}h) = {vial_activity_at_compounding:.2f} mCi")
    print(f"   - Current concentration = {vial_activity_at_compounding:.2f} mCi / {VIAL_VOLUME_mL:.1f} mL = {vial_concentration_at_compounding:.2f} mCi/mL\n")

    print("4. Final Volume to be drawn into the syringe:")
    print("   - Formula: Volume = Activity Needed / Current Concentration")
    # Final equation with numbers
    print(f"   - Volume = {activity_needed_at_compounding:.2f} mCi / {vial_concentration_at_compounding:.2f} mCi/mL")
    print("\n---------------------------------------------------------------------")
    print(f"The volume of Indium 111 chloride that must be drawn is: {volume_to_draw_ml:.2f} mL")
    print("---------------------------------------------------------------------")


solve_radiopharmacy_calculation()