def identify_animal():
  """
  Identifies the animal in the image and prints its scientific name.
  """
  is_animal_present = True
  animal_common_name = "longhorn beetle"
  animal_scientific_name = "Stictoleptura cordigera"

  if is_animal_present:
    print(f"Yes, there is an animal in the image. It is a type of {animal_common_name}.")
    print(f"The animal's scientific name is: {animal_scientific_name}")
  else:
    print("No animal was identified in the image.")

identify_animal()