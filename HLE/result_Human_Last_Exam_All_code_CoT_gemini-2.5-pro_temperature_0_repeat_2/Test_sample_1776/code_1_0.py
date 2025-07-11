import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    """
    # --- Step 1: Define constants and initial values ---
    A_inject = 8.0  # mCi, desired activity at injection
    V_cal = 10.0  # mL, stock calibration volume
    A_cal = 10.0  # mCi, stock calibration activity
    T_half_days = 2.8 # days, half-life of Indium-111
    T_half_hours = T_half_days * 24 # 67.2 hours

    # Define the relevant times using the provided dates
    injection_time = datetime(2024, 12, 23, 8, 0)
    compounding_time = datetime(2024, 12, 23, 4, 4)
    calibration_time = datetime(2024, 12, 26, 12, 0)

    # --- Step 2: Calculate the decay constant (lambda) ---
    # The formula is: lambda = ln(2) / T_half
    decay_constant_hourly = math.log(2) / T_half_hours

    # --- Step 3: Calculate the activity needed at compounding time ---
    # This accounts for decay between compounding and injection.
    # The formula is: A_needed = A_final * e^(lambda * t)
    
    # Time delta between compounding and injection
    t_decay_dose_delta = injection_time - compounding_time
    t_decay_dose_hours = t_decay_dose_delta.total_seconds() / 3600
    
    # Activity needed at compounding time
    A_compound = A_inject * math.exp(decay_constant_hourly * t_decay_dose_hours)

    # --- Step 4: Calculate the activity and concentration of the stock vial at compounding time ---
    # This is a reverse decay calculation from the calibration time.
    # The formula is: A_initial = A_final * e^(lambda * t)

    # Time delta between compounding and calibration
    t_pre_decay_stock_delta = calibration_time - compounding_time
    t_pre_decay_stock_hours = t_pre_decay_stock_delta.total_seconds() / 3600

    # Activity of the stock vial at compounding time
    A_stock_now = A_cal * math.exp(decay_constant_hourly * t_pre_decay_stock_hours)
    
    # Concentration of the stock vial at compounding time
    C_stock_now = A_stock_now / V_cal

    # --- Step 5: Calculate the final volume to draw ---
    # The formula is: Volume = Activity_needed / Concentration
    V_draw = A_compound / C_stock_now

    # --- Print the detailed explanation and results ---
    print("--- Radiopharmacy Calculation ---")
    print(f"\n1. First, we determine the activity needed at the compounding time (4:04 am).")
    print(f"   - The dose of {A_inject:.1f} mCi is required at injection (8:00 am).")
    print(f"   - The time between compounding and injection is {t_decay_dose_hours:.2f} hours.")
    print(f"   - To account for decay, the activity needed at compounding is:")
    print(f"     A_needed = {A_inject:.1f} mCi * e^({decay_constant_hourly:.6f} * {t_decay_dose_hours:.2f}h) = {A_compound:.3f} mCi")

    print(f"\n2. Next, we determine the concentration of the stock vial at compounding time.")
    print(f"   - The vial is calibrated for {A_cal:.1f} mCi in {V_cal:.1f} mL at 12:00 pm on Dec 26.")
    print(f"   - The time from compounding (4:04 am, Dec 23) to calibration is {t_pre_decay_stock_hours:.2f} hours.")
    print(f"   - The vial's activity at compounding time was:")
    print(f"     A_stock = {A_cal:.1f} mCi * e^({decay_constant_hourly:.6f} * {t_pre_decay_stock_hours:.2f}h) = {A_stock_now:.3f} mCi")
    print(f"   - The vial's concentration at compounding time was:")
    print(f"     C_stock = {A_stock_now:.3f} mCi / {V_cal:.1f} mL = {C_stock_now:.3f} mCi/mL")

    print(f"\n3. Finally, we calculate the volume to draw.")
    print(f"   - The required volume is the needed activity divided by the stock concentration.")
    print("\n--- Final Equation ---")
    print(f"   Volume = {A_compound:.3f} mCi / {C_stock_now:.3f} mCi/mL")
    print(f"\n   The volume to draw is {V_draw:.2f} mL.")
    
    return V_draw

if __name__ == '__main__':
    solve_radiopharmacy_calculation()