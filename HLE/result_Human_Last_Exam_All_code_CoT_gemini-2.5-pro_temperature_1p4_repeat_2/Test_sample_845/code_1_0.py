def find_hrv_term():
    """
    This function identifies and explains the term used for the time
    value between heartbeats in HRV analysis.
    """
    # The value between heartbeats is called the Inter-Beat Interval.
    term = "Inter-Beat Interval (IBI)"
    
    # Another common term, specific to ECG, is R-R Interval.
    # Since the context is general (PPG), IBI is the more appropriate term.
    
    # Formulate the answer as a clear sentence.
    explanation = f"The value computed between consecutive heartbeats for HRV analysis is called the {term}. It represents the time duration between each beat, typically measured in milliseconds."
    
    print(explanation)

find_hrv_term()