def get_hrv_value_name():
    """
    This function defines and prints the name of the value measured between
    heartbeats for HRV analysis.
    """
    # The value is the time interval between successive heartbeats.
    # The most common technical term is "R-R interval".
    # Another widely used term, especially with PPG sensors, is "Inter-Beat Interval (IBI)".
    value_name = "R-R interval"
    alternate_name = "Inter-Beat Interval (IBI)"

    print(f"The value computed between heartbeats to track HRV is called the '{value_name}' (or '{alternate_name}').")

# Execute the function to print the answer.
get_hrv_value_name()