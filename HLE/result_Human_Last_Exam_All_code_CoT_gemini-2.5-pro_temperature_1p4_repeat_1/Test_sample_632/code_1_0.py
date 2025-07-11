def find_removed_feature_info(year):
  """
  This function sets up the parameters for the historical query.
  """
  location = "Piazza della Rotonda"
  event_description = "A unique architectural feature was removed."
  
  # The problem specifies the year of removal.
  # As per the instructions, we will print the number involved in the final output.
  print(f"Investigating a historical event in {location}.")
  print(f"Description: {event_description}")
  print(f"The removal occurred in the year: {year}")

# Execute the function with the year from the question.
find_removed_feature_info(1950)