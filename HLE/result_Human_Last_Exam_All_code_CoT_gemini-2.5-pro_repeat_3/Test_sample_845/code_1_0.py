def get_hrv_value_name():
  """
  This function provides the name of the value computed between heartbeats for HRV.
  """
  # The value is the time interval between consecutive heartbeats.
  # The most common technical term is "R-R interval".
  value_name = "R-R interval"
  print(f"The value computed between heartbeats for HRV is called the '{value_name}'.")

get_hrv_value_name()