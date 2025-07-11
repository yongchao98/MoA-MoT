import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose.
    The code follows these steps:
    1. Define all known constants and times.
    2. Calculate the time difference between compounding and injection (t1).
    3. Calculate the activity needed at compounding time (A_required) to decay to the final dose.
    4. Calculate the time difference between compounding and calibration (t2).
    5. Calculate the total activity in the stock vial at compounding time (A_vial_now).
    6. Calculate the concentration of the stock vial at compounding time (C_vial_now).
    7. Calculate the final volume to draw by dividing required activity by the vial concentration.
    """

    # --- Step 1: Define Constants and Times ---
    # Half-life of Indium-111 in hours
    T_HALF_LIFE_DAYS = 2.805
    T_HALF_LIFE_HOURS = T_HALF_LIFE_DAYS * 24

    # Decay constant (lambda) in hours^-1
    LAMBDA_H = math.log(2) / T_HALF_LIFE_HOURS

    # Dose and vial information
    A_dose = 8.0  # mCi, activity required for injection
    A_cal = 10.0  # mCi, activity of the vial at calibration time
    V_cal = 10.0  # mL, total volume of the vial

    # Relevant times
    time_compound = datetime(2024, 12, 23, 4, 4)
    time_injection = datetime(2024, 12, 23, 8, 0)
    time_calibration = datetime(2024, 12, 26, 12, 0)

    # --- Step 2: Calculate Activity Needed at Compounding Time (`A_required`) ---
    # This is the activity needed at 4:04 am to have 8 mCi at 8:00 am.
    # We calculate the time difference (t1) between compounding and injection.
    t1_delta = time_injection - time_compound
    t1_hours = t1_delta.total_seconds() / 3600

    # A_required = A_dose * e^(lambda * t1)
    A_required = A_dose * math.exp(LAMBDA_H * t1_hours)
    
    print("--- Step 1: Calculate Activity Required at Compounding Time ---")
    print(f"The half-life of In-111 is {T_HALF_LIFE_DAYS} days or {T_HALF_LIFE_HOURS:.2f} hours.")
    print(f"The decay constant (λ) is ln(2) / {T_HALF_LIFE_HOURS:.2f} h = {LAMBDA_H:.6f} h⁻¹.")
    print(f"Time from compounding (4:04 am) to injection (8:00 am) is {t1_hours:.2f} hours.")
    print(f"Required Activity = {A_dose} mCi * e^({LAMBDA_H:.6f} * {t1_hours:.2f})")
    print(f"Required Activity = {A_required:.3f} mCi\n")


    # --- Step 3: Calculate Vial Concentration at Compounding Time (`C_vial_now`) ---
    # We calculate the vial's activity at 4:04 am based on its calibration at a later time.
    # We calculate the time difference (t2) from compounding back to calibration.
    t2_delta = time_calibration - time_compound
    t2_hours = t2_delta.total_seconds() / 3600
    
    # A_vial_now = A_cal * e^(lambda * t2)
    A_vial_now = A_cal * math.exp(LAMBDA_H * t2_hours)
    
    # C_vial_now = A_vial_now / V_cal
    C_vial_now = A_vial_now / V_cal
    
    print("--- Step 2: Calculate Vial Concentration at Compounding Time ---")
    print(f"Time from compounding (Dec 23, 4:04 am) to calibration (Dec 26, 12:00 pm) is {t2_hours:.2f} hours.")
    print(f"Total vial activity at compounding = {A_cal} mCi * e^({LAMBDA_H:.6f} * {t2_hours:.2f})")
    print(f"Total vial activity at compounding = {A_vial_now:.3f} mCi")
    print(f"Vial concentration = {A_vial_now:.3f} mCi / {V_cal} mL = {C_vial_now:.3f} mCi/mL\n")

    # --- Step 4: Calculate Final Volume to Draw ---
    V_draw = A_required / C_vial_now
    
    print("--- Step 3: Calculate Final Volume to Draw ---")
    print("Volume to Draw = Required Activity / Vial Concentration")
    print(f"Volume to Draw = {A_required:.3f} mCi / {C_vial_now:.3f} mCi/mL")
    print(f"\nThe volume of Indium-111 chloride that must be drawn is: {V_draw:.3f} mL")
    
    return V_draw

# Execute the function and capture the final answer for submission format
final_volume = solve_radiopharmacy_calculation()
# Final answer format as requested
print(f"<<<{final_volume:.3f}>>>")
