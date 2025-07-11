def calculate_map():
    """
    This function calculates the Mean Arterial Pressure (MAP) based on the patient's resting blood pressure.
    """
    # Patient's resting vital signs from the case
    systolic_bp = 135  # mmHg
    diastolic_bp = 85   # mmHg

    # The formula for MAP is: Diastolic BP + 1/3 * (Systolic BP - Diastolic BP)
    pulse_pressure = systolic_bp - diastolic_bp
    mean_arterial_pressure = diastolic_bp + (1/3) * pulse_pressure

    print("To fulfill the request for a calculation, we can calculate the patient's Mean Arterial Pressure (MAP) at rest.")
    print("The formula is: MAP = Diastolic BP + 1/3 * (Systolic BP - Diastolic BP)\n")
    print("Using the patient's values:")
    print(f"MAP = {diastolic_bp} + 1/3 * ({systolic_bp} - {diastolic_bp})")
    print(f"MAP = {diastolic_bp} + 1/3 * ({pulse_pressure})")
    print(f"The patient's resting MAP is approximately {mean_arterial_pressure:.2f} mmHg.")

calculate_map()