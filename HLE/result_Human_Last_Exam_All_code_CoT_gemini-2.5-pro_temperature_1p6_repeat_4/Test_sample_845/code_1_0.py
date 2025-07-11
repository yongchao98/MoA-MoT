def get_hrv_value_name():
    """
    This function explains and returns the name of the value computed between heartbeats for HRV.
    """
    # HRV is the measure of the variation in time between consecutive heartbeats.
    # This time duration is known as the R-R interval (the interval between two consecutive R-waves of an ECG)
    # or more generally as the inter-beat interval (IBI).
    value_name = "R-R interval"
    print(f"The value computed as a time duration between heartbeats to calculate HRV is called the: {value_name}")

if __name__ == "__main__":
    get_hrv_value_name()