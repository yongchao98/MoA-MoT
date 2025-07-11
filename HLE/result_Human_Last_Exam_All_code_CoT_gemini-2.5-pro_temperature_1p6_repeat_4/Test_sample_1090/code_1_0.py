def solve_clinical_case():
    """
    Analyzes the patient's case to determine the most appropriate next step in management.
    """
    
    # Patient Data and Key Calculations
    resting_hr = 76
    standing_hr = 118
    hr_increase = standing_hr - resting_hr
    
    full_knee_extension_angle = 180
    limited_knee_extension_angle = 165
    extension_deficit = full_knee_extension_angle - limited_knee_extension_angle

    print("Step 1: Analyzing the patient's key clinical findings.")
    print(f"The patient's heart rate increases by {standing_hr} - {resting_hr} = {hr_increase} bpm upon standing, indicating significant orthostatic intolerance and deconditioning.")
    print(f"The key physical finding is a knee extension deficit of {full_knee_extension_angle} - {limited_knee_extension_angle} = {extension_deficit} degrees due to resistance.")
    print("-" * 20)

    print("Step 2: Identifying the primary barrier to recovery.")
    print("The patient was ambulatory with a cane prior to admission but now cannot walk despite physical therapy.")
    print("The resistance to passive knee extension is a classic sign of post-stroke spasticity, which has likely worsened due to his recent severe illness and immobility.")
    print("This spasticity is the primary mechanical barrier preventing effective participation in physical therapy and a return to ambulation.")
    print("-" * 20)
    
    print("Step 3: Determining the most appropriate management step.")
    print("While deconditioning, orthostatic intolerance, and poor nutritional status are contributing factors, addressing the spasticity is the most critical next step to restore function.")
    print("A Physical Medicine and Rehabilitation (PM&R) specialist (physiatrist) is an expert in managing complex rehabilitation needs, including post-stroke spasticity.")
    print("A PM&R consultation will lead to a comprehensive plan that could include specialized therapy techniques, bracing, and/or medical management of spasticity.")
    print("-" * 20)

    print("Conclusion: The single most appropriate next step is to involve the specialist best equipped to manage this problem.")
    print("Final Answer: Consult Physical Medicine and Rehabilitation.")

solve_clinical_case()
<<<Consult Physical Medicine and Rehabilitation>>>