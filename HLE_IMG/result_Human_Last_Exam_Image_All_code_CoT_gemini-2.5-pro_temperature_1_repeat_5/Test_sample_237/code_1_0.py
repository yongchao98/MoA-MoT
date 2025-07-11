def identify_animal():
  """
  This function identifies the animal in the image and prints its scientific name.
  """
  is_animal_present = True
  animal_type = "Longhorn Beetle (Insect)"
  scientific_name = "Stenocorus meridianus"

  if is_animal_present:
    print("Yes, there is an animal in the image.")
    print(f"The animal is a: {animal_type}")
    print(f"Its scientific name is: {scientific_name}")
  else:
    print("No animal was identified in the image.")

identify_animal()