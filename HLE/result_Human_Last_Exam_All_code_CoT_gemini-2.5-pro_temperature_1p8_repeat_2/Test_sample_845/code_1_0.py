def get_hrv_value_name():
    """
    This function explains and returns the name of the value used to compute HRV.
    HRV (Heart Rate Variability) is the measure of the variation in time between
    consecutive heartbeats. This time interval itself is the fundamental value.
    """
    # The most common term for the time interval between consecutive heartbeats,
    # originating from ECG (electrocardiogram) readings, is the "R-R interval".
    # It represents the time between two successive R-waves.
    answer = "R-R interval"

    print("Heart Rate Variability (HRV) is calculated from the time intervals between successive heartbeats.")
    print(f"This beat-to-beat interval is most commonly called the: {answer}")
    print("(Another closely related term, especially for PPG-based sensors, is 'Inter-Beat Interval' or IBI).")

get_hrv_value_name()