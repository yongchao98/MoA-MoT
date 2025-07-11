def get_hrv_value_name():
    """
    This function explains and prints the name of the value used for HRV computation.
    """
    # HRV is the variation in the time between consecutive heartbeats.
    # The value itself is the measurement of this time interval.
    # While PPG (Photoplethysmography) measures pulse-to-pulse intervals, the terminology
    # is often borrowed from ECG (Electrocardiogram) terminology.
    value_name = "R-R interval"
    
    # In ECG, the most prominent spike is the 'R' wave. The time between two
    # consecutive R waves is the "R-R interval". This is the most common term.
    # A more general term is "Inter-beat Interval" or IBI.
    
    print(f"The value computed between heartbeats, which is the time interval between them, is called the '{value_name}'.")
    print("A more general term, often used with PPG, is the 'Inter-beat Interval' (IBI).")

if __name__ == "__main__":
    get_hrv_value_name()