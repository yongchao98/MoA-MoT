def identify_animal():
  """
  Identifies the animal in the image and provides its scientific name.
  """
  animal_present = True
  common_name = "Red-brown Longhorn Beetle"
  scientific_name = "Stictoleptura rubra"

  if animal_present:
    print(f"Yes, there is an animal in the image.")
    print(f"It is commonly known as the {common_name}.")
    print(f"The animal's scientific name is: {scientific_name}")
  else:
    print("No animal was identified in the image.")

identify_animal()