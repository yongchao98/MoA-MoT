import sys

def identify_organism():
  """
  Identifies the black organism in the image.
  """
  organism_name = "Spinellus fusiger (a type of pin mold)"
  description = (
      "The small black organisms visible on the surface of the larger shelf fungus are a "
      "mycoparasitic fungus, meaning a fungus that grows on another fungus. "
      "Their appearance, with tiny black-headed pins, is characteristic of a pin mold. "
      "The most likely identification is Spinellus fusiger, also known as the Bonnet Mould."
  )
  
  print(f"Name: {organism_name}")
  print(f"Description: {description}")

identify_organism()