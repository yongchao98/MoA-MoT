def get_hrv_value_name():
    """
    This function explains and returns the name of the value measured between
    heartbeats for HRV calculation.
    """
    # The value captured between heartbeats is the time interval.
    # While PPG measures Pulse-to-Pulse Interval (PPI), the most standard
    # term used in HRV analysis, derived from ECG terminology, is the R-R interval.
    # It represents the time between the peaks of two consecutive R waves in an ECG signal.
    # Inter-beat Interval (IBI) is another common and more general term.

    answer = "R-R interval"
    explanation = "The value between heartbeats used to compute HRV is the time duration between them. The most common term for this is the R-R interval (the time between two successive R-waves of the QRS signal on an ECG), or more generally, the Inter-beat Interval (IBI)."
    
    print(explanation)
    print(f"\nThe specific value is called: {answer}")

if __name__ == "__main__":
    get_hrv_value_name()