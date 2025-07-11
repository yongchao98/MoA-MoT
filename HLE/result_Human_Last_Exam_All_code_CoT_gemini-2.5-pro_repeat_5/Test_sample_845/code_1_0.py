def get_hrv_interval_name():
  """
  This function provides the name of the fundamental value used to compute HRV.
  HRV is calculated from the variations in the time between successive heartbeats.
  """
  # The term for the time interval between consecutive heartbeats.
  # This term originates from ECG readings (R-peak to R-peak) but is
  # the standard term used for beat-to-beat intervals from any source, including PPG.
  interval_name = "R-R interval"
  
  print(f"The value representing the time between heartbeats, used to compute HRV, is called the: {interval_name}")

get_hrv_interval_name()