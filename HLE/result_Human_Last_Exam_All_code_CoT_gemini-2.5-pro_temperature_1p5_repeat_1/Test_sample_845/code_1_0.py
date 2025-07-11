def find_hrv_value_name():
    """
    Explains and identifies the name of the value computed between heartbeats for HRV.
    """
    print("Heart Rate Variability (HRV) is calculated based on the precise time intervals between successive heartbeats.")
    print("This beat-to-beat interval is the fundamental measurement.")
    print("In clinical and scientific literature, especially when derived from an ECG (Electrocardiogram), this interval is named after the 'R' wave, which is the most prominent peak in a heartbeat's electrical signal.")
    print("\nThe value computed between heartbeats is called the:")
    
    # The term is composed of two letters and a word.
    term_part_1 = "R"
    term_part_2 = "R"
    term_part_3 = "interval"
    
    # We will print the full term now.
    print(f"{term_part_1}-{term_part_2} {term_part_3}")
    
    print("\nNote: When using PPG (Photoplethysmography), a more generic term is 'Inter-Beat Interval' (IBI), but 'R-R interval' is the most commonly used term in the context of HRV.")

if __name__ == '__main__':
    find_hrv_value_name()