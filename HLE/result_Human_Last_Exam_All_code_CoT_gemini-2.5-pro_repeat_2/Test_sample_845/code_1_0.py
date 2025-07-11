def print_hrv_term():
  """
  This function defines and prints the term for the value computed between heartbeats for HRV analysis.
  """
  # The value between heartbeats is the time interval.
  term = "Inter-Beat Interval (IBI)"
  
  # Print the explanation.
  print(f"The value computed as the time between heartbeats to measure HRV is called the: {term}")

# Execute the function to display the answer.
print_hrv_term()