def analyze_patient_vitals():
    """
    This function calculates the change in the patient's heart rate from resting to standing
    and prints the details of the calculation.
    """
    # Patient's heart rate data from the clinical case
    resting_hr = 76  # beats/min
    standing_hr = 118 # beats/min

    # Calculate the change in heart rate
    hr_change = standing_hr - resting_hr

    # Print the findings and the calculation
    print("A key finding is the significant tachycardia upon standing.")
    print("The calculation for the change in heart rate is as follows:")
    print(f"{standing_hr} (standing) - {resting_hr} (resting) = {hr_change} beats/min")
    print(f"This increase of {hr_change} beats/min highlights the patient's significant deconditioning and physiologic stress.")

analyze_patient_vitals()