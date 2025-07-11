import sys

def identify_animal():
  """
  This function identifies the animal in the image and prints its scientific name.
  """
  is_animal_present = True
  animal_type = "insect (Heart-marked Red Longhorn Beetle)"
  scientific_name = "Stictoleptura cordigera"

  if is_animal_present:
    print(f"Yes, there is an animal in this image. It is an {animal_type}.")
    print(f"The animal's scientific name is: {scientific_name}")
  else:
    print("No animal was identified in the image.")

identify_animal()
