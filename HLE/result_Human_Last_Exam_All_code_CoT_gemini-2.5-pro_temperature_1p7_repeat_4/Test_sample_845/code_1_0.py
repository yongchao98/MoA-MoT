def get_hrv_value_name():
  """
  Explains and prints the name of the value used to compute HRV.
  HRV is the variation in the time between heartbeats.
  This time is most commonly called the R-R interval.
  """
  # The name of the value is the "R-R interval".
  # It represents the time between the peaks of two consecutive 'R' waves on an ECG.
  # When using PPG, it is sometimes called the Inter-Beat Interval (IBI), but R-R interval is the standard term in HRV literature.
  answer = "The value computed between heartbeats for HRV is called the R-R interval or Inter-Beat Interval (IBI)."
  print(answer)

get_hrv_value_name()