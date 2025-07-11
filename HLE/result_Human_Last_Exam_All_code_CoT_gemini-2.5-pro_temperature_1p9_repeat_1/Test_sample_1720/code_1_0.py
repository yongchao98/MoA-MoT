def analyze_vitals(systolic_bp, diastolic_bp):
    """
    Calculates Mean Arterial Pressure (MAP) to assess tissue perfusion.
    """
    # The formula for Mean Arterial Pressure (MAP) is: DBP + 1/3 * (SBP - DBP)
    # A MAP < 65 mmHg generally indicates inadequate organ perfusion.
    pulse_pressure = systolic_bp - diastolic_bp
    mean_arterial_pressure = diastolic_bp + (1/3) * pulse_pressure

    print("Analyzing patient's hemodynamic status based on blood pressure.")
    print(f"Given Systolic BP: {systolic_bp}")
    print(f"Given Diastolic BP: {diastolic_bp}")
    print("\nCalculating Mean Arterial Pressure (MAP):")
    print(f"Equation: MAP = Diastolic BP + (1/3) * (Systolic BP - Diastolic BP)")
    print(f"Step 1: MAP = {diastolic_bp} + (1/3) * ({systolic_bp} - {diastolic_bp})")
    print(f"Step 2: MAP = {diastolic_bp} + (1/3) * ({pulse_pressure})")
    print(f"Step 3: MAP = {diastolic_bp} + {round((1/3) * pulse_pressure, 2)}")
    print(f"Final Result: MAP = {round(mean_arterial_pressure, 2)} mmHg")
    if mean_arterial_pressure < 65:
        print("\nConclusion: The patient's MAP is dangerously low, confirming a state of shock requiring immediate, aggressive treatment.")
    else:
        print("\nConclusion: The patient's MAP is borderline low, indicating hypotension that requires immediate, aggressive treatment to prevent worsening shock.")


# Vital signs from the patient case
patient_systolic = 90
patient_diastolic = 60

analyze_vitals(patient_systolic, patient_diastolic)