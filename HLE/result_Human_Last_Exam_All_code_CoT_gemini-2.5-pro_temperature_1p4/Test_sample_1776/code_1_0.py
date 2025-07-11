import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium 111 chloride for a patient dose.
    """
    # --- Given Information ---
    
    # Radionuclide properties
    half_life_days = 2.805  # Half-life of Indium-111 in days

    # Desired dose information
    desired_activity = 8.0  # mCi at time of injection
    injection_time_str = "2024-12-23 08:00:00"

    # Compounding information (when the dose is drawn)
    compounding_time_str = "2024-12-23 04:04:00"

    # Stock vial information
    stock_cal_activity = 10.0  # mCi at calibration time
    stock_cal_volume = 10.0    # mL
    stock_cal_time_str = "2024-12-26 12:00:00"

    # --- Step-by-Step Calculation ---

    print("Step-by-step calculation for the required volume of Indium 111 chloride:")
    print("-" * 70)

    # 1. Calculate the decay constant (lambda)
    lambda_val = math.log(2) / half_life_days
    print(f"1. The decay constant (λ) for Indium-111 is calculated as:")
    print(f"   λ = ln(2) / {half_life_days:.3f} days = {lambda_val:.5f} days⁻¹\n")

    # Convert string times to datetime objects for accurate calculations
    injection_time = datetime.strptime(injection_time_str, "%Y-%m-%d %H:%M:%S")
    compounding_time = datetime.strptime(compounding_time_str, "%Y-%m-%d %H:%M:%S")
    stock_cal_time = datetime.strptime(stock_cal_time_str, "%Y-%m-%d %H:%M:%S")

    # 2. Calculate the required activity at the time of compounding
    # This accounts for the decay between compounding and injection.
    time_delta_decay = injection_time - compounding_time
    time_delta_decay_days = time_delta_decay.total_seconds() / (24 * 3600)
    activity_at_compounding = desired_activity * math.exp(lambda_val * time_delta_decay_days)
    
    print(f"2. Calculate the activity required at compounding time ({compounding_time.strftime('%I:%M %p, %b %d')}):")
    print(f"   Time between compounding and injection = {time_delta_decay_days:.4f} days")
    print(f"   Required Activity = {desired_activity:.2f} mCi * e^({lambda_val:.5f} * {time_delta_decay_days:.4f})")
    print(f"   Required Activity = {activity_at_compounding:.4f} mCi\n")

    # 3. Calculate the concentration of the stock vial at compounding time
    # This accounts for decay from the future calibration time back to the compounding time.
    time_delta_stock = stock_cal_time - compounding_time
    time_delta_stock_days = time_delta_stock.total_seconds() / (24 * 3600)
    stock_activity_at_compounding = stock_cal_activity * math.exp(lambda_val * time_delta_stock_days)
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_cal_volume

    print(f"3. Calculate the stock vial concentration at compounding time:")
    print(f"   Time from compounding to vial calibration = {time_delta_stock_days:.4f} days")
    print(f"   Stock Activity = {stock_cal_activity:.2f} mCi * e^({lambda_val:.5f} * {time_delta_stock_days:.4f}) = {stock_activity_at_compounding:.4f} mCi")
    print(f"   Stock Concentration = {stock_activity_at_compounding:.4f} mCi / {stock_cal_volume:.1f} mL = {stock_concentration_at_compounding:.4f} mCi/mL\n")

    # 4. Calculate the final volume to draw
    volume_to_draw = activity_at_compounding / stock_concentration_at_compounding

    print("4. Final volume calculation:")
    print("   Volume to Draw = Required Activity / Stock Concentration")
    # This prints the final equation with the calculated numbers
    print(f"   Volume to Draw = {activity_at_compounding:.4f} mCi / {stock_concentration_at_compounding:.4f} mCi/mL")
    print("-" * 70)
    print(f"\nThe volume of Indium 111 chloride that must be drawn is: {volume_to_draw:.2f} mL")
    
    # Return final answer for capture
    return f"<<<{volume_to_draw:.2f}>>>"

# Run the calculation and print the final result in the specified format
final_answer = solve_radiopharmacy_calculation()
print(final_answer)
