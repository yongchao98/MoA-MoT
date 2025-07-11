def get_hrv_value_name():
    """
    This function explains and provides the name for the value measured between
    heartbeats to calculate Heart Rate Variability (HRV).
    """
    # The value is the time elapsed between two successive heartbeats.
    # While technically called Pulse-to-Pulse Interval (PPI) when using PPG,
    # the most common and universally understood term is Inter-beat Interval (IBI).
    # R-R interval is also frequently used, originating from ECG measurements.
    value_name = "Inter-beat interval (IBI)"
    
    print(f"The value representing the time between heartbeats, which is used to compute HRV, is most commonly called the '{value_name}'.")

# Execute the function to print the answer.
get_hrv_value_name()