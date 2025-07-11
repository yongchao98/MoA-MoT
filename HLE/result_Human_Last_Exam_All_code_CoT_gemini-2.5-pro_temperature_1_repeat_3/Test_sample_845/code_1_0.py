def get_hrv_term():
  """
  This function explains and returns the name for the value measured between heartbeats.
  """
  term_name = "Inter-Beat Interval (IBI)"
  alternative_name = "R-R Interval"
  explanation = (
      f"The value computed between heartbeats for HRV analysis is called the '{term_name}'.\n"
      f"When measured using an ECG (electrocardiogram), it is specifically known as the '{alternative_name}'."
  )
  return explanation

# Print the explanation
print(get_hrv_term())