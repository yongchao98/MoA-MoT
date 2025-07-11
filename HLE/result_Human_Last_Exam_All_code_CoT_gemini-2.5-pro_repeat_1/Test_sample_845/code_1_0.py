def get_hrv_value_name():
    """
    This function defines and prints the name of the value computed between heartbeats for HRV.
    """
    # The value is the time measured between consecutive heartbeats.
    value_name = "Inter-Beat Interval"
    abbreviation = "IBI"

    # Print the explanation and the name of the value.
    print(f"HRV (Heart Rate Variability) is computed from the variation in the time between consecutive heartbeats.")
    print(f"The value representing the time between each beat is called the: {value_name} ({abbreviation}).")
    print("When using an ECG, this is often specifically called the R-R interval.")

# Execute the function to display the answer.
get_hrv_value_name()