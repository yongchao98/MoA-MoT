import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium 111 chloride for a patient dose.
    """
    # --- Constants and Initial Values ---
    half_life_days = 2.80  # Half-life of Indium-111 in days
    vial_activity_cal = 10.0  # mCi at calibration
    vial_volume_cal = 10.0  # mL
    dose_activity_injection = 8.0  # mCi at injection time

    # --- Timestamps ---
    time_calibration = datetime(2024, 12, 26, 12, 0)
    time_compounding = datetime(2024, 12, 23, 4, 4)
    time_injection = datetime(2024, 12, 23, 8, 0)

    # --- Calculations ---
    # Calculate the decay constant (lambda) in units of per day
    lambda_per_day = math.log(2) / half_life_days

    # --- Step 1: Calculate the activity needed at compounding time (4:04 am) ---
    print("--- Step 1: Calculate Dose Activity Needed at Compounding Time ---")
    
    # Time from compounding to injection
    decay_time_to_injection = time_injection - time_compounding
    decay_time_to_injection_days = decay_time_to_injection.total_seconds() / (24 * 3600)
    
    # We need to find the initial activity (A₀) that will decay to 8.0 mCi.
    # Formula: A₀ = A * e^(λ * t)
    activity_needed_at_compounding = dose_activity_injection * math.exp(lambda_per_day * decay_time_to_injection_days)

    print(f"The dose must be {dose_activity_injection:.2f} mCi at injection time ({time_injection.strftime('%I:%M %p')}).")
    print(f"We need to find the required activity at compounding time ({time_compounding.strftime('%I:%M %p')}).")
    print(f"Time for decay between compounding and injection: {decay_time_to_injection} ({decay_time_to_injection_days:.4f} days).")
    print(f"The decay equation is: Activity_at_compounding = Activity_at_injection * e^(λ * t)")
    print("The final equation with numbers plugged in is:")
    print(f"Activity_at_compounding = {dose_activity_injection:.2f} mCi * e^({lambda_per_day:.4f} * {decay_time_to_injection_days:.4f})")
    print(f"Result: The activity that must be drawn at 4:04 am is {activity_needed_at_compounding:.4f} mCi.\n")

    # --- Step 2: Calculate the vial's concentration at compounding time (4:04 am) ---
    print("--- Step 2: Calculate Vial Concentration at Compounding Time ---")
    
    # Time from compounding (past) to calibration (future)
    time_from_comp_to_cal = time_calibration - time_compounding
    time_from_comp_to_cal_days = time_from_comp_to_cal.total_seconds() / (24 * 3600)

    # Calculate the vial's activity at the earlier compounding time.
    # Formula: A₀ = A * e^(λ * t)
    vial_activity_at_compounding = vial_activity_cal * math.exp(lambda_per_day * time_from_comp_to_cal_days)
    vial_concentration_at_compounding = vial_activity_at_compounding / vial_volume_cal

    print(f"The vial is calibrated for {vial_activity_cal:.2f} mCi at {time_calibration.strftime('%Y-%m-%d %I:%M %p')}.")
    print(f"Time from compounding to calibration: {time_from_comp_to_cal} ({time_from_comp_to_cal_days:.4f} days).")
    print("To find the activity at the earlier compounding time, the decay equation is: Activity_at_compounding = Activity_at_calibration * e^(λ * t)")
    print("The final equation with numbers plugged in is:")
    print(f"Vial_Activity_at_compounding = {vial_activity_cal:.2f} mCi * e^({lambda_per_day:.4f} * {time_from_comp_to_cal_days:.4f})")
    print(f"Result: The vial's total activity at 4:04 am is {vial_activity_at_compounding:.4f} mCi.")
    print(f"The vial's concentration is {vial_activity_at_compounding:.4f} mCi / {vial_volume_cal:.2f} mL = {vial_concentration_at_compounding:.4f} mCi/mL.\n")

    # --- Step 3: Calculate the final volume to draw ---
    print("--- Step 3: Calculate Final Volume to Draw ---")
    volume_to_draw = activity_needed_at_compounding / vial_concentration_at_compounding

    print("To find the required volume, we divide the required activity (Step 1) by the vial's concentration (Step 2).")
    print("The final equation with numbers plugged in is:")
    print(f"Volume = {activity_needed_at_compounding:.4f} mCi / {vial_concentration_at_compounding:.4f} mCi/mL")
    print("-" * 20)
    print(f"The final volume to be drawn is: {volume_to_draw:.2f} mL")
    print("-" * 20)
    
    return volume_to_draw

if __name__ == "__main__":
    final_volume = solve_radiopharmacy_calculation()
    # The final answer is wrapped in <<<>>> as requested.
    # Using round(value, 2) to present a practical answer for a syringe.
    print(f"\n<<< {round(final_volume, 2)} >>>")
