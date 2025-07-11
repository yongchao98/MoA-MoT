def explain_hrv_value():
    """
    Explains the term for the value measured between heartbeats for HRV analysis.
    """
    explanation = (
        "The value computed between heartbeats to measure Heart Rate Variability (HRV) is the time duration between those consecutive beats.\n\n"
        "This value is most commonly called the 'Inter-beat Interval' (IBI).\n\n"
        "In the context of an ECG (electrocardiogram), where heartbeats are identified by 'R' waves, this is also famously known as the 'R-R interval'."
    )
    print(explanation)

if __name__ == "__main__":
    explain_hrv_value()