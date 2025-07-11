def find_hrv_value_name():
    """
    This function explains and prints the name of the value computed between
    heartbeats to measure Heart Rate Variability (HRV).
    """
    # The value captured between consecutive heartbeats is the time interval.
    # This is the fundamental data point for all HRV analysis.
    value_name = "Inter-Beat Interval (IBI)"
    alternative_name = "R-R Interval" # Specifically from ECG QRS complex

    # The prompt mentions PPG, for which IBI is a very common term.
    # R-R interval is also widely used, even when the source isn't strictly an ECG.
    
    print(f"The value computed between heartbeats for HRV is the time interval between them.")
    print(f"This is most commonly called the '{value_name}' or '{alternative_name}'.")

find_hrv_value_name()