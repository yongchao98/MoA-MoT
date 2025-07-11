def identify_animal_in_image():
  """
  This function provides the identification of the animal in the image.
  """
  animal_present = True
  common_name = "Ruddy longhorn beetle"
  scientific_name = "Anastrangalia sanguinolenta"

  if animal_present:
    print(f"Yes, there is an animal in the image. It is a {common_name}.")
    print(f"Its scientific name is: {scientific_name}")
  else:
    print("No animal could be identified in the image.")

identify_animal_in_image()