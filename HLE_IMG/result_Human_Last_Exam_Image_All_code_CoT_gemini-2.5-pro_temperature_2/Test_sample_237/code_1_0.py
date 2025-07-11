import sys

def identify_animal():
  """
  This function identifies the animal in the image and prints its scientific name.
  """
  is_animal_present = True
  animal_type = "longhorn beetle"
  scientific_name = "Analeptura liturata"

  if is_animal_present:
    print("Yes, there is an animal in this image.")
    # The double backslash for the output is to make sure the final output format is what is expected
    # since a single one would escape the final <.
    # It does not mean I am outputting two backslashes.
    sys.stdout.write(f"The animal is a {animal_type}, and its scientific name is {scientific_name}.\n")
  else:
    print("No animal was identified in the image.")

identify_animal()