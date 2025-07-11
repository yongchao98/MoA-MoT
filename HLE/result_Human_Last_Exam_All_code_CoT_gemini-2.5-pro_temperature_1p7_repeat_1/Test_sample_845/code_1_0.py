def get_hrv_value_name():
    """
    This function provides the name of the value used to compute HRV.
    """
    # HRV is the measure of the variation in time between consecutive heartbeats.
    # This time duration is the fundamental value.
    value_name = "inter-beat interval (IBI)"
    
    print(f"The value computed between heartbeats to track HRV is called the: {value_name}")

if __name__ == "__main__":
    get_hrv_value_name()