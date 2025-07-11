def find_hrv_value_name():
    """
    This function identifies and explains the term for the value computed
    between heartbeats for HRV analysis.
    """
    # The term for the interval between successive R-peaks (heartbeats)
    # is the "R-R interval". It's the fundamental unit for calculating HRV.
    term = "R-R interval"
    
    # Another common name for this value is the "inter-beat interval" or IBI.
    alternative_term = "inter-beat interval (IBI)"
    
    print(f"The value representing the time between consecutive heartbeats, used to compute HRV, is called the '{term}'.")
    print(f"It is also commonly referred to as the '{alternative_term}'.")

# Execute the function to provide the answer
find_hrv_value_name()