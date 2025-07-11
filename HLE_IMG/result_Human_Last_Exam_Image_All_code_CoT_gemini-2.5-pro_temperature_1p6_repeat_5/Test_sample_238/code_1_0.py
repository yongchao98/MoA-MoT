import sys

def identify_organism():
  """
  This function identifies the organisms in the provided image, focusing on the black one.
  """
  primary_fungus = "Bracket fungus (likely a Trametes species)"
  secondary_organism_description = "Tiny, black-tipped, pin-like structures"
  secondary_organism_name = "Pin mold (a type of mycoparasitic fungus)"

  print(f"The image contains several organisms.")
  print(f"The large, shelf-like organism is a {primary_fungus}.")
  print(f"The {secondary_organism_description} growing on the larger fungus are a {secondary_organism_name}.")
  print(f"\nThe black organism in the photo is the: {secondary_organism_name}")

if __name__ == "__main__":
  identify_organism()
