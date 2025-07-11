def classify_microscope_image():
  """
  Provides a classification for the microorganism shown in the image.
  """
  morphology = "The microorganism displays dimorphism, with two forms visible: round-to-oval budding yeast cells (blastoconidia) and chains of elongated cells forming pseudohyphae (distinguished by constrictions at the septa)."
  organism_type = "This morphology is characteristic of a fungus, specifically a yeast belonging to the genus Candida. The presence of both yeasts and pseudohyphae is a hallmark of an invasive infection, often by Candida albicans."
  stain = "The stain appears to be Periodic acid-Schiff (PAS), which stains the fungal cell walls a distinct magenta color."

  print("Classification based on the image:")
  print(f"\n1. Morphology: {morphology}")
  print(f"\n2. Type of Microorganism: {organism_type}")
  print(f"\n3. Probable Stain: {stain}")

classify_microscope_image()