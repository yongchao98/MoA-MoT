def get_hrv_term():
    """
    This function defines and prints the name of the value used to compute HRV.
    """
    # The value is the time interval between consecutive heartbeats.
    # There are two common names for this interval.
    primary_term = "R-R interval"
    general_term = "Inter-Beat Interval (IBI)"

    print(f"The value measured between heartbeats to compute HRV is the time interval itself.")
    print(f"This is most commonly called the '{primary_term}' or more generally the '{general_term}'.")

get_hrv_term()