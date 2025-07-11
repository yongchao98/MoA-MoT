def find_atlantic_indicator_species():
  """
  This function identifies the species in the photo that indicates
  an Atlantic rather than a Pacific coast location.

  The most identifiable species with a limited geographic range in the photo
  is the common periwinkle.

  Scientific Name: Littorina littorea
  Distribution: Native to the Northeast Atlantic (Europe). Introduced and
                now widespread on the North American Atlantic coast. It is
                not established on the Pacific coast.

  Therefore, its presence strongly suggests the photo was taken on the
  Atlantic coast.
  """
  scientific_name = "Littorina littorea"
  print(f"The scientific name of the indicator species is: {scientific_name}")

find_atlantic_indicator_species()