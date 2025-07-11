import sys

def get_bud_for_single_dose_ampule():
  """
  Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule
  opened in a sterile environment according to USP <797> guidelines.
  """
  
  # Define the parameters based on the user's query
  container_type = "single-dose container (ampule or vial)"
  environment = "ISO Class 5 or better sterile environment"
  
  # According to USP <797> (revised 2023), the BUD for a single-dose container
  # once opened or punctured in an ISO 5 environment is 12 hours.
  bud_in_hours = 12
  
  # Print the explanation and the final answer
  print(f"The question is about the Beyond-Use Date (BUD) for a {container_type}.")
  print(f"The container is opened (punctured or broken) in a {environment}.")
  print("Based on USP <797> guidelines, the BUD is determined by these factors.")
  print("\nThe final answer is:")
  print(f"The BUD is {bud_in_hours} hours from the time of puncture.")

# Execute the function
if __name__ == '__main__':
  get_bud_for_single_dose_ampule()