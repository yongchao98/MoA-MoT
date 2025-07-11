def find_hrv_value_name():
    """
    This function explains and identifies the name of the value
    measured between heartbeats for HRV calculation.
    """

    # Description of the concepts
    hrv_context = """
HRV (Heart Rate Variability) is a measure of the variation in time between successive heartbeats.
This variation is not random; it reflects the regulation of the heart by the autonomic nervous system.
To calculate HRV, we first need to measure the precise time between each beat.
"""
    print(hrv_context)

    # The name of the value
    value_name_explanation = "The value measured between heartbeats is called the 'Inter-beat interval' or IBI. When measured using an electrocardiogram (ECG), it is more specifically known as the 'R-R interval', named after the 'R' peak of the QRS complex in the ECG waveform."
    print(value_name_explanation)

    final_answer_term = "Inter-beat interval (IBI)"
    print(f"\nTherefore, the value is called: {final_answer_term}")


if __name__ == "__main__":
    find_hrv_value_name()