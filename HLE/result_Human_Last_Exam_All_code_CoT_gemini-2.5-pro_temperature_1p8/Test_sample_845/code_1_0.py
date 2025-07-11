def get_hrv_value_name():
    """
    This function explains the name of the value used to compute HRV.
    HRV (Heart Rate Variability) is the measure of the variation in time
    between consecutive heartbeats. This time duration is the fundamental
    value used in all HRV calculations.
    """
    value_name_general = "inter-beat interval (IBI)"
    value_name_ecg = "R-R interval"

    print(f"The value computed between heartbeats for HRV analysis is the time duration between them.")
    print(f"This is most commonly called the '{value_name_general}'.")
    print(f"When measured with an ECG, it is specifically referred to as the '{value_name_ecg}'.")

if __name__ == "__main__":
    get_hrv_value_name()