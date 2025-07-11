def get_hrv_term():
    """
    This function explains the term used for the value computed between heartbeats for HRV.
    """
    # HRV is the measure of the variation in time between consecutive heartbeats.
    # This specific time value measured between beats is most commonly called the R-R interval.
    # Even when using PPG, where it is technically the Inter-Beat Interval (IBI),
    # the term R-R interval is widely used.
    
    answer = "The value computed between heartbeats for HRV analysis is called the R-R interval (or Inter-Beat Interval - IBI)."
    print(answer)

get_hrv_term()