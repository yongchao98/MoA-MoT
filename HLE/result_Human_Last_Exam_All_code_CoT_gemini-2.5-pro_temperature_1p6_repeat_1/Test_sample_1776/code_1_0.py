import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 for a patient dose.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in days
    T_half_days = 2.805
    # Required dose activity at injection time in mCi
    activity_at_injection = 8.0
    # Stock vial calibration activity in mCi
    stock_cal_activity = 10.0
    # Stock vial total volume in mL
    stock_cal_volume = 10.0

    # --- Timestamps ---
    # Compounding time: Dec 23, 4:04 am
    compounding_time = datetime(2024, 12, 23, 4, 4)
    # Injection time: Dec 23, 8:00 am
    injection_time = datetime(2024, 12, 23, 8, 0)
    # Calibration time: Dec 26, 12:00 pm
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # --- Calculations ---
    # Convert half-life to hours and calculate the decay constant (lambda)
    half_life_hours = T_half_days * 24.0
    decay_constant = math.log(2) / half_life_hours

    # Step 1: Calculate activity needed at compounding time
    # Time elapsed from compounding to injection
    dt_comp_to_inj = (injection_time - compounding_time).total_seconds() / 3600.0  # in hours
    # Use pre-decay formula: A_initial = A_final * e^(lambda * t)
    activity_needed_at_compounding = activity_at_injection * math.exp(decay_constant * dt_comp_to_inj)

    # Step 2: Calculate the stock vial's concentration at compounding time
    # Time elapsed from compounding to calibration (stock activity decays over this time)
    dt_comp_to_cal = (calibration_time - compounding_time).total_seconds() / 3600.0  # in hours
    # Use pre-decay formula to find activity at the earlier time (compounding)
    stock_activity_at_compounding = stock_cal_activity * math.exp(decay_constant * dt_comp_to_cal)
    # Calculate concentration
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_cal_volume

    # Step 3: Calculate the volume to draw
    # Volume = Activity_Needed / Concentration
    volume_to_draw = activity_needed_at_compounding / stock_concentration_at_compounding

    # --- Print Results with Equations ---
    print("This script calculates the volume of Indium-111 chloride to be drawn.\n")

    print("Step 1: Calculate the activity needed at compounding time (4:04 am).")
    print("This accounts for decay between compounding and injection.")
    print(f"A_needed = (Activity at injection) * e^(lambda * t)")
    print(f"A_needed = {activity_at_injection:.2f} mCi * e^({decay_constant:.6f} * {dt_comp_to_inj:.2f} hours) = {activity_needed_at_compounding:.2f} mCi\n")

    print("Step 2: Calculate the stock vial concentration at compounding time.")
    print("This accounts for the decay from the time of compounding to the future calibration time.")
    print(f"A_stock = (Activity at calibration) * e^(lambda * t)")
    print(f"A_stock_at_compounding = {stock_cal_activity:.2f} mCi * e^({decay_constant:.6f} * {dt_comp_to_cal:.2f} hours) = {stock_activity_at_compounding:.2f} mCi")
    print(f"C_stock = (Activity_stock) / (Volume_stock)")
    print(f"C_stock_at_compounding = {stock_activity_at_compounding:.2f} mCi / {stock_cal_volume:.2f} mL = {stock_concentration_at_compounding:.2f} mCi/mL\n")

    print("Step 3: Calculate the final volume to draw.")
    print("Volume = (Activity Needed) / (Stock Concentration)")
    print(f"Volume to Draw = {activity_needed_at_compounding:.2f} mCi / {stock_concentration_at_compounding:.2f} mCi/mL = {volume_to_draw:.2f} mL")


solve_radiopharmacy_calculation()