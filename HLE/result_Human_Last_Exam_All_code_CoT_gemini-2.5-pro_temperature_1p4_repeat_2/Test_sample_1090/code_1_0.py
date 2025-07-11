def calculate_hr_change():
    """
    This function calculates the change in heart rate from resting to standing
    based on the values provided in the clinical case.
    """
    
    # Vital signs from the case text
    resting_hr = 76  # beats/min at rest
    standing_hr = 118 # beats/min when attempting to stand

    # Calculate the increase in heart rate
    hr_increase = standing_hr - resting_hr

    # Print the equation as requested
    print(f"The patient's heart rate increases significantly upon attempting to stand. The calculation is:")
    print(f"{standing_hr} - {resting_hr} = {hr_increase}")

calculate_hr_change()