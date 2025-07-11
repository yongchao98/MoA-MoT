def identify_animal():
  """
  This function identifies the animal in the image and prints its scientific name.
  """
  animal_present = True
  common_name = "Longhorn Beetle (likely a Heart-bearing longhorn beetle)"
  scientific_name = "Stictoleptura cordigera"

  if animal_present:
    print("Yes, there is an animal in the image.")
    print(f"It appears to be a type of {common_name}.")
    print(f"The scientific name is: {scientific_name}")
  else:
    print("No animal could be identified in the image.")

identify_animal()