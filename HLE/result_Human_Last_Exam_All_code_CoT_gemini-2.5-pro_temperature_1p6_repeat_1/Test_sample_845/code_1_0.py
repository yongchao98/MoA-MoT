def find_hrv_value_name():
    """
    This function explains and prints the name of the value computed between heartbeats for HRV.
    """
    # HRV is the measure of the variation in time between consecutive heartbeats.
    # The value itself is the duration of this time interval.
    value_name = "R-R interval"
    explanation = (
        "The value computed between heartbeats for calculating HRV is the time interval between them. "
        f"This is most commonly called the '{value_name}'."
    )
    print(explanation)
    print("\nThe term 'R-R interval' originates from ECG (electrocardiogram) readings, where it represents the time between two successive R-peaks.")
    print("When using PPG (photoplethysmography), the equivalent term is Inter-Beat Interval (IBI), but 'R-R interval' is often used interchangeably.")

find_hrv_value_name()