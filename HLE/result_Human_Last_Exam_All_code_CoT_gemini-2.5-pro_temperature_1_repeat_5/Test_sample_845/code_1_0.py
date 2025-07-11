def get_hrv_value_name():
  """
  Explains the value used to compute Heart Rate Variability (HRV).
  """
  # HRV is computed from the variation in the time between consecutive heartbeats.
  value_name = "R-R interval"
  alternative_name = "inter-beat interval (IBI)"
  
  print(f"The value computed between heartbeats for HRV is the time duration between them.")
  print(f"This is most commonly called the '{value_name}' or '{alternative_name}'.")

get_hrv_value_name()