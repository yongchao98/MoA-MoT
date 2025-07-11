def get_hrv_term():
  """
  Provides the name for the value computed between heartbeats for HRV.
  """
  # The value is the time interval between consecutive heartbeats.
  # This is most commonly called the R-R interval, a term derived from ECGs,
  # but is also known as the Inter-Beat Interval (IBI).
  answer = "The value computed between heartbeats is the time interval, which is called the R-R interval or Inter-Beat Interval (IBI)."
  print(answer)

get_hrv_term()