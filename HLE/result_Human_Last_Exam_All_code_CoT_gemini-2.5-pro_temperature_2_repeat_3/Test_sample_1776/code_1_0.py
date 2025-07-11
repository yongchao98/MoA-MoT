import math
from datetime import datetime, timedelta

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium 111 chloride needed for a patient dose.
    """
    # --- Step 1: Define Constants and Times ---
    T_half_days = 2.8
    T_half_hours = T_half_days * 24

    A_inject_target = 8.0  # mCi, desired activity at injection time
    A_cal = 10.0  # mCi, vial's calibrated activity
    V_cal = 10.0  # mL, vial's calibrated volume

    t_cal = datetime(2024, 12, 26, 12, 0)
    t_compound = datetime(2024, 12, 23, 4, 4)
    t_inject = datetime(2024, 12, 23, 8, 0)

    # --- Step 2: Calculate Decay Constant (λ) ---
    lambda_val = math.log(2) / T_half_hours

    # --- Step 3: Calculate Time Differences in Hours ---
    # Time from compounding to future calibration
    delta_t_vial = t_cal - t_compound
    delta_t_vial_hours = delta_t_vial.total_seconds() / 3600

    # Time from compounding to injection
    delta_t_dose = t_inject - t_compound
    delta_t_dose_hours = delta_t_dose.total_seconds() / 3600
    
    print("Plan:")
    print("1. Calculate the dose activity needed at the time of compounding.")
    print("2. Calculate the vial's concentration at the time of compounding.")
    print("3. Divide the required dose activity by the vial's concentration to find the volume.\n")

    # --- Step 4: Calculate Activity Needed for the Dose at Compounding Time ---
    # We need to draw a dose that will decay to A_inject_target at t_inject.
    # A_draw = A_inject_target * e^(λ * Δt_dose)
    A_draw = A_inject_target * math.exp(lambda_val * delta_t_dose_hours)
    
    print("--- Calculation for the Dose to be Drawn ---")
    print(f"The dose must be {A_inject_target} mCi at injection time ({t_inject.strftime('%I:%M %p')}).")
    print(f"The time between compounding ({t_compound.strftime('%I:%M %p')}) and injection is {delta_t_dose_hours:.2f} hours.")
    print(f"Equation: Activity_to_draw = {A_inject_target:.1f} mCi * e^({lambda_val:.4f} * {delta_t_dose_hours:.2f} hours)")
    print(f"Required activity at compounding time: {A_draw:.4f} mCi\n")


    # --- Step 5: Calculate Vial's Activity Concentration at Compounding Time ---
    # We need the vial's activity at t_compound, based on its future calibration.
    # A_vial_at_compound = A_cal * e^(λ * Δt_vial)
    A_vial_at_compound = A_cal * math.exp(lambda_val * delta_t_vial_hours)
    C_vial_at_compound = A_vial_at_compound / V_cal
    
    print("--- Calculation for the Vial's Concentration ---")
    print(f"The vial is calibrated for {A_cal:.1f} mCi at {t_cal.strftime('%Y-%m-%d %I:%M %p')}.")
    print(f"The time between compounding and calibration is {delta_t_vial_hours:.2f} hours.")
    print(f"Equation: Vial_Activity = {A_cal:.1f} mCi * e^({lambda_val:.4f} * {delta_t_vial_hours:.2f} hours)")
    print(f"Vial's total activity at compounding time: {A_vial_at_compound:.4f} mCi")
    print(f"Vial's concentration = {A_vial_at_compound:.4f} mCi / {V_cal:.1f} mL = {C_vial_at_compound:.4f} mCi/mL\n")

    # --- Step 6: Calculate Final Volume ---
    # Volume = Activity needed / Concentration
    V_draw = A_draw / C_vial_at_compound
    
    print("--- Final Volume Calculation ---")
    print(f"Equation: Volume_to_draw = Required_Activity / Vial_Concentration")
    print(f"Volume to draw = {A_draw:.4f} mCi / {C_vial_at_compound:.4f} mCi/mL")
    print(f"The final volume to draw into the syringe is: {V_draw:.2f} mL")
    
    return V_draw

if __name__ == "__main__":
    final_volume = solve_radiopharmacy_calculation()
    # The final answer is wrapped for the system.
    # print(f"\n<<<final_answer>>>\n{final_volume:.2f}\n<<<final_answer>>>")

#<<<3.65>>>