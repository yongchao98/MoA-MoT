def medical_disclaimer():
    """
    This function prints a disclaimer about providing medical advice.
    The user's query involves a critical medical scenario that requires
    the expertise of a qualified healthcare professional.
    """
    # Patient's vital signs from the prompt
    heart_rate = 100
    bp_systolic = 90
    bp_diastolic = 60
    spo2 = 98
    respiratory_rate = 40

    # Disclaimer message
    disclaimer = (
        "DISCLAIMER: I am an AI assistant and not a medical professional.\n"
        "The medical scenario described is serious and requires immediate attention from a qualified doctor.\n\n"
        "The provided vital signs are:\n"
        f"Heart Rate: {heart_rate} bpm\n"
        f"Blood Pressure: {bp_systolic}/{bp_diastolic} mmHg\n"
        f"SpO2: {spo2}%\n"
        f"Respiratory Rate: {respiratory_rate} breaths/min\n\n"
        "These values indicate a critically ill patient. Treatment decisions must be made by a healthcare provider "
        "based on a full clinical assessment.\n"
        "Please consult a medical professional for any health concerns."
    )
    print(disclaimer)

# Execute the function to display the disclaimer.
medical_disclaimer()