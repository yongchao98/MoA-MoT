def explain_hrv_term():
  """
  Explains and identifies the term for the beat-to-beat interval used in HRV calculations.
  """
  explanation = "Heart Rate Variability (HRV) measures the variation in time between successive heartbeats."
  primary_term = "R-R interval"
  context = "This term originates from ECG (Electrocardiogram) measurements, where the interval is measured between two consecutive R-peaks."
  secondary_term = "Inter-Beat Interval (IBI)"
  ppg_context = "In the context of PPG (Photoplethysmography), this value is also precisely called the Inter-Beat Interval (IBI). However, 'R-R interval' is the most widely used term in HRV literature regardless of the measurement technology."

  print(explanation)
  print(f"The value, or time interval, computed between heartbeats is called the '{primary_term}'.")
  print(f"{context}")
  print(f"It is also known as the '{secondary_term}'.")
  print(f"{ppg_context}")

explain_hrv_term()