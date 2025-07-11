def get_hrv_value_name():
    """
    This function provides the name of the value used to compute HRV.
    """
    # HRV is the measure of variation in time between successive heartbeats.
    # This time is the core value used for all HRV calculations.
    value_name_1 = "R-R interval"
    value_name_2 = "Inter-Beat Interval (IBI)"

    # Both terms are correct and often used interchangeably.
    # R-R interval comes from ECG terminology (the 'R' peak).
    # IBI is a more general term.
    print(f"The value is called the {value_name_1} or {value_name_2}.")

if __name__ == "__main__":
    get_hrv_value_name()