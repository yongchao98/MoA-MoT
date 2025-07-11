def get_hrv_interval_name():
  """
  Provides the standard name for the interval between heartbeats used in HRV calculation.
  """
  # The value computed between heartbeats for HRV analysis has a specific name.
  # This value is the time duration between two consecutive heartbeats.
  # In ECG, this is the time between two 'R' waves, the most prominent peaks.
  # This term is widely adopted even when using other methods like PPG.
  hrv_value_name = "R-R interval"

  print(f"The value computed between heartbeats to calculate HRV is called the: {hrv_value_name}")

get_hrv_interval_name()