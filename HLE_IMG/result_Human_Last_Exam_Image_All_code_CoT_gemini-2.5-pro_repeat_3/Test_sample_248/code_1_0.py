def find_atlantic_indicator_species():
  """
  This function identifies and prints the scientific name of the species in the image
  that is found on the Atlantic coast but not the Pacific coast.
  """
  # The image contains several species, including barnacles and periwinkles.
  # The Common Periwinkle (Littorina littorea) is native to the European Atlantic
  # and introduced to the North American Atlantic. It is not found on the Pacific coast.
  # Its presence is therefore a definitive indicator of an Atlantic location.
  scientific_name = "Littorina littorea"
  print(f"The scientific name of the species indicating an Atlantic coast location is: {scientific_name}")

find_atlantic_indicator_species()