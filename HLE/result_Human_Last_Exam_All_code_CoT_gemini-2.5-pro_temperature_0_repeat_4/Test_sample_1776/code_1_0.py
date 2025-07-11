import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Problem Data ---
    # Dose ordered for the patient
    activity_needed_at_injection = 8.0  # mCi

    # Time of injection
    injection_time = datetime(2024, 12, 23, 8, 0)

    # Time of kit compounding
    compounding_time = datetime(2024, 12, 23, 4, 4)

    # Stock vial information
    stock_calibration_activity = 10.0  # mCi
    stock_calibration_volume = 10.0    # mL
    stock_calibration_time = datetime(2024, 12, 26, 12, 0)

    # Physical constant for Indium-111
    half_life_days = 2.8
    half_life_hours = half_life_days * 24

    # --- Step-by-step Calculation ---

    print("This script calculates the required volume of Indium-111 chloride.\n")

    # Step 1: Calculate the decay constant (lambda) for Indium-111
    # lambda = ln(2) / T_half
    decay_constant_h = math.log(2) / half_life_hours
    print(f"Step 1: Calculate the decay constant (λ) for Indium-111")
    print(f"Half-life (T½) = {half_life_days} days = {half_life_hours:.1f} hours")
    print(f"Decay constant (λ) = ln(2) / {half_life_hours:.1f} h = {decay_constant_h:.6f} per hour\n")

    # Step 2: Calculate the activity needed at the time of compounding.
    # This accounts for the decay between compounding and injection.
    # A_compound = A_injection * e^(λ * t1)
    time_delta_compounding_to_injection = injection_time - compounding_time
    t1_hours = time_delta_compounding_to_injection.total_seconds() / 3600
    activity_needed_at_compounding = activity_needed_at_injection * math.exp(decay_constant_h * t1_hours)
    print(f"Step 2: Calculate the activity needed at compounding time ({compounding_time.strftime('%I:%M %p %b %d')})")
    print(f"The dose of {activity_needed_at_injection:.1f} mCi is needed at {injection_time.strftime('%I:%M %p %b %d')}.")
    print(f"Time between compounding and injection = {t1_hours:.3f} hours")
    print(f"Activity needed at compounding = {activity_needed_at_injection:.1f} mCi * e^({decay_constant_h:.6f} * {t1_hours:.3f}) = {activity_needed_at_compounding:.3f} mCi\n")

    # Step 3: Calculate the activity of the stock vial at the time of compounding.
    # This accounts for the decay from the future calibration date back to the compounding time.
    # A_stock_now = A_calibration * e^(λ * t2)
    time_delta_compounding_to_calibration = stock_calibration_time - compounding_time
    t2_hours = time_delta_compounding_to_calibration.total_seconds() / 3600
    stock_activity_at_compounding = stock_calibration_activity * math.exp(decay_constant_h * t2_hours)
    print(f"Step 3: Calculate the stock vial's activity at compounding time")
    print(f"The stock vial is calibrated for {stock_calibration_activity:.1f} mCi at {stock_calibration_time.strftime('%I:%M %p %b %d')}.")
    print(f"Time from compounding to calibration = {t2_hours:.3f} hours")
    print(f"Stock activity at compounding = {stock_calibration_activity:.1f} mCi * e^({decay_constant_h:.6f} * {t2_hours:.3f}) = {stock_activity_at_compounding:.3f} mCi\n")

    # Step 4: Calculate the concentration of the stock vial at the time of compounding.
    # Concentration = Activity / Volume
    stock_concentration_at_compounding = stock_activity_at_compounding / stock_calibration_volume
    print(f"Step 4: Calculate the stock vial's concentration at compounding time")
    print(f"Stock concentration = {stock_activity_at_compounding:.3f} mCi / {stock_calibration_volume:.1f} mL = {stock_concentration_at_compounding:.3f} mCi/mL\n")

    # Step 5: Calculate the volume of Indium-111 chloride to draw.
    # Volume = Desired Activity / Concentration
    volume_to_draw = activity_needed_at_compounding / stock_concentration_at_compounding
    print(f"Step 5: Calculate the final volume to draw")
    print(f"The final equation is: Volume = (Activity needed at compounding) / (Concentration at compounding)")
    print(f"Volume = {activity_needed_at_compounding:.3f} mCi / {stock_concentration_at_compounding:.3f} mCi/mL")
    print(f"\n-----------------------------------------------------------------")
    print(f"The required volume to draw at {compounding_time.strftime('%I:%M %p')} is: {volume_to_draw:.2f} mL")
    print(f"-----------------------------------------------------------------")

if __name__ == '__main__':
    solve_radiopharmacy_calculation()