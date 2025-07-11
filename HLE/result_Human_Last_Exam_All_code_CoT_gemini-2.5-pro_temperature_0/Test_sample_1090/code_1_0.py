def analyze_patient_status():
    """
    This function analyzes key patient metrics to determine the most likely cause
    of his failure to progress in physical therapy.
    """
    # Patient's vital signs
    resting_heart_rate = 76  # bpm
    standing_heart_rate = 118 # bpm

    # Patient's Body Mass Index
    bmi = 18.5  # kg/m2

    # Calculate the increase in heart rate upon standing
    heart_rate_increase = standing_heart_rate - resting_heart_rate

    print("Analyzing the patient's failure to ambulate:")
    print(f"The patient's heart rate at rest is {resting_heart_rate} bpm.")
    print(f"When attempting to stand, his heart rate increases to {standing_heart_rate} bpm.")
    print(f"The calculation for the heart rate increase is: {standing_heart_rate} - {resting_heart_rate} = {heart_rate_increase} bpm.")
    print(f"This large increase in heart rate with minimal exertion, combined with a low BMI of {bmi} kg/m2, strongly suggests severe deconditioning due to a catabolic state and malnutrition.")
    print("\nConclusion: The patient's nutritional status must be addressed before he can be expected to succeed in physical therapy.")

analyze_patient_status()