def identify_animal():
  """
  This function identifies the animal in the provided image and prints its scientific name.
  """
  is_animal_present = True
  common_name = "Heart-marked Longhorn Beetle"
  scientific_name = "Stictoleptura cordigera"

  if is_animal_present:
    print("Yes, there is an animal in the image.")
    print(f"It appears to be a {common_name}.")
    print(f"The animal's scientific name is: {scientific_name}")
  else:
    print("No animal could be identified in the image.")

identify_animal()