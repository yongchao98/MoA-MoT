def classify_microorganism_from_image():
  """
  This function prints a detailed classification of the microorganism
  based on the provided microscopic image.
  """
  morphology = "The image shows two distinct morphological forms: ovoid, budding yeast cells (blastoconidia) and filamentous structures known as pseudohyphae. The pseudohyphae are formed from chains of elongated yeast cells that remain attached after budding and are characterized by constrictions at the septations."
  microorganism_type = "The presence of both budding yeasts and pseudohyphae is a hallmark characteristic of a fungus, specifically a yeast belonging to the genus *Candida*. This presentation is typical of a *Candida* infection (candidiasis) in tissue."
  stain_info = "The staining method is likely Periodic acid-Schiff (PAS), which stains the polysaccharide-rich fungal cell walls a bright magenta/purple, contrasting with the blue/green counterstained host tissue."

  print("Classification of the Observed Microorganism:")
  print("-" * 40)
  print(f"Morphology: {morphology}")
  print(f"\nType of Microorganism: {microorganism_type}")
  print(f"\nProbable Stain: {stain_info}")
  print("-" * 40)

classify_microorganism_from_image()