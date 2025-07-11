def get_hrv_value_name():
    """
    This function explains and prints the name of the value used to compute HRV.
    """
    # HRV is the variation in the time between consecutive heartbeats.
    # This time measurement is the fundamental value used for HRV analysis.
    value_name = "R-R interval"
    alternative_name = "Inter-Beat Interval (IBI)"

    print(f"The value computed between heartbeats to calculate Heart Rate Variability (HRV) is called the '{value_name}'.")
    print(f"It is also commonly referred to as the '{alternative_name}'.")

# Execute the function to provide the answer.
get_hrv_value_name()