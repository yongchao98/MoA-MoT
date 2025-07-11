def get_hrv_term():
    """
    This function provides the name for the value computed between heartbeats to measure HRV.
    """
    # The value computed between heartbeats is the time interval.
    # When measured using PPG, this is commonly called the Inter-beat interval (IBI).
    # When measured using an electrocardiogram (ECG), it is called the R-R interval.
    term_name = "Inter-beat interval (IBI)"
    alternative_term = "R-R interval"
    
    print(f"The value computed between heartbeats to calculate HRV is called the '{term_name}'.")
    print(f"It is the time duration between consecutive heartbeats. Another common name, especially in the context of ECGs, is the '{alternative_term}'.")

get_hrv_term()