def identify_animal():
  """
  This function identifies the animal in the image and prints its scientific name.
  """
  animal_present = True
  animal_type = "insect (Longhorn Beetle)"
  scientific_name = "Stictoleptura cordigera"

  if animal_present:
    print("Yes, there is an animal in this image.")
    print(f"The animal is an {animal_type}.")
    print(f"What is the animal's scientific name? \nIts scientific name is: {scientific_name}")

identify_animal()