def get_hrv_value_name():
  """
  Provides the name for the value computed between heartbeats for HRV.
  """
  # The value is the time duration between two consecutive heartbeats.
  # For PPG, this is commonly called the Inter-Beat Interval (IBI).
  # For ECG, it is called the R-R Interval.
  answer = "The value computed between heartbeats for HRV analysis is called the Inter-Beat Interval (IBI) or R-R Interval."
  print(answer)

get_hrv_value_name()