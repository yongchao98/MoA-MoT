def identify_beetle():
  """
  This script provides the identification of the animal in the image.
  """
  animal_present = True
  common_name = "Black-and-red longhorn beetle"
  scientific_name = "Stenurella melanura"

  if animal_present:
    print("Yes, there is an animal in the image.")
    print(f"The animal is an insect, likely a {common_name}.")
    print(f"Its scientific name is probably {scientific_name}.")
  else:
    print("No animal identified in the image.")

identify_beetle()