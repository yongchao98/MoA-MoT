import math
from datetime import datetime, timedelta

def calculate_volume_to_draw():
    """
    Calculates the volume of Indium-111 Chloride required for a patient dose.
    """
    # Step 1: Define all known parameters from the problem
    desired_dose_activity = 8.0  # mCi
    vial_cal_activity = 10.0  # mCi
    vial_cal_volume = 10.0  # mL
    
    # Half-life of Indium-111 in hours
    t_half_hours = 2.805 * 24  # 67.32 hours

    # Define the relevant date and time points
    # The year 2024 is used as specified. Timezone is assumed to be consistent.
    cal_time = datetime(2024, 12, 26, 12, 0)  # 12:00 PM on Dec 26
    inject_time = datetime(2024, 12, 23, 8, 0)  # 8:00 AM on Dec 23
    compound_time = datetime(2024, 12, 23, 4, 4) # 4:04 AM on Dec 23 (for context)

    print("Problem Analysis:")
    print(f"Desired dose: {desired_dose_activity} mCi at {inject_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"Vial calibration: {vial_cal_activity} mCi in {vial_cal_volume} mL at {cal_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"Compounding time: {compound_time.strftime('%Y-%m-%d %I:%M %p')}")
    print(f"Half-life of Indium-111: {t_half_hours:.2f} hours")
    print("-" * 30)

    # Step 2: Calculate the decay constant (lambda)
    lambda_val = math.log(2) / t_half_hours
    print("Step 1: Calculate the decay constant (lambda) for Indium-111.")
    print(f"lambda = ln(2) / T_half_life")
    print(f"lambda = {math.log(2):.6f} / {t_half_hours:.2f} hours = {lambda_val:.6f} per hour\n")
    
    # Step 3: Calculate the time difference between injection and calibration
    time_delta = inject_time - cal_time
    time_delta_hours = time_delta.total_seconds() / 3600
    print("Step 2: Calculate time decay factor from calibration time to injection time.")
    print(f"Time difference = Injection Time - Calibration Time")
    print(f"Time difference = {inject_time} - {cal_time} = {time_delta_hours:.2f} hours\n")

    # Step 4: Calculate the initial concentration of the vial at calibration time
    initial_concentration = vial_cal_activity / vial_cal_volume
    print("Step 3: Calculate the initial concentration of the vial at calibration time.")
    print(f"Initial Concentration = Vial Activity / Vial Volume")
    print(f"Initial Concentration = {vial_cal_activity} mCi / {vial_cal_volume} mL = {initial_concentration:.2f} mCi/mL\n")

    # Step 5: Calculate the required volume
    # Formula: Volume = (Desired_Dose / Initial_Concentration) * e^(lambda * time_delta_hours)
    # This calculates the volume of stock solution needed at calibration time which will decay to the right dose concentration
    # and since we are drawing from the same stock vial, this volume remains constant.
    required_volume = (desired_dose_activity / initial_concentration) * math.exp(lambda_val * time_delta_hours)
    
    print("Step 4: Calculate the volume to draw using the radioactive decay formula.")
    print("The volume to draw is calculated based on what volume of the vial (at calibration) would decay to the correct dose at the injection time.")
    print("Formula: Volume = (Dose / Initial_Conc) * e^(lambda * time_difference_hours)")
    print("\nFinal Equation:")
    print(f"Volume (mL) = ({desired_dose_activity} mCi / {initial_concentration:.1f} mCi/mL) * e^({lambda_val:.6f} * {time_delta_hours:.2f})")
    print(f"Volume (mL) = {desired_dose_activity / initial_concentration:.2f} * e^({lambda_val * time_delta_hours:.4f})")
    print(f"Volume (mL) = {desired_dose_activity / initial_concentration:.2f} * {math.exp(lambda_val * time_delta_hours):.4f}")
    
    print("\nResult:")
    print(f"The required volume to draw into the syringe is {required_volume:.2f} mL.")

    return required_volume

if __name__ == '__main__':
    final_volume = calculate_volume_to_draw()
    # The final answer is wrapped according to the instructions
    # print(f'<<<{final_volume:.2f}>>>')