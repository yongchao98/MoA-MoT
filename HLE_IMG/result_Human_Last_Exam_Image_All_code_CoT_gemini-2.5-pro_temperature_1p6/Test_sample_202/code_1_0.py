def identify_rock_specimens():
  """
  This function identifies the rock specimens based on their visual characteristics.
  H: Reddish, vesicular (full of holes) - typical of Scoria.
  I: Composed of angular rock fragments cemented together - a sedimentary Breccia.
  J: Light-colored, highly vesicular, low-density volcanic rock - Pumice.
  K: Very coarse-grained with large, interlocking crystals - an igneous Pegmatite.
  """
  
  identifications = {
      "H": "Scoria",
      "I": "Breccia",
      "J": "Pumice",
      "K": "Pegmatite"
  }
  
  print("Identification of Rock Specimens:")
  for specimen_letter, rock_name in identifications.items():
    print(f"Specimen {specimen_letter} is {rock_name}.")

identify_rock_specimens()