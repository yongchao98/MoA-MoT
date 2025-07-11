def analyze_patient_data():
    """
    This script analyzes key numerical data from the patient's case
    to determine the most likely cause of his inability to ambulate and
    suggest the next step in management.
    """

    # Resting vitals
    resting_hr = 76
    resting_sbp = 135
    resting_dbp = 85

    # Vitals upon standing/walking attempt
    standing_hr = 118
    standing_sbp = 150
    standing_dbp = 90

    # Neurological finding
    full_knee_extension_angle = 180
    patient_knee_extension_angle = 165

    # Calculations
    hr_change = standing_hr - resting_hr
    sbp_change = standing_sbp - resting_sbp
    dbp_change = standing_dbp - resting_dbp
    knee_extension_deficit = full_knee_extension_angle - patient_knee_extension_angle

    print("Patient Clinical Data Analysis:")
    print("-" * 35)
    print("Change in Vitals upon Standing:")
    print(f"Heart Rate Increase: {standing_hr} - {resting_hr} = {hr_change} bpm")
    print(f"Systolic BP Change: {standing_sbp} - {resting_sbp} = {sbp_change} mmHg")
    print(f"Diastolic BP Change: {standing_dbp} - {resting_dbp} = {dbp_change} mmHg")
    print("\nKey Physical Exam Finding:")
    print(f"Knee Extension Deficit: {full_knee_extension_angle}° - {patient_knee_extension_angle}° = {knee_extension_deficit}°")
    print("-" * 35)
    
    print("\nClinical Conclusion:")
    print("The patient's inability to ambulate is best explained by a combination of hospital-acquired deconditioning and a new left knee flexion contracture, identified by the 15° extension deficit on his paretic side. Orthostatic hypotension is ruled out.")
    
    print("\nMost Appropriate Next Step:")
    print("The primary issue is a new musculoskeletal impairment. The management should therefore focus on rehabilitation.")
    
# Execute the analysis
analyze_patient_data()

# Final answer is provided below after the code block.