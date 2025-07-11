def get_species_indices():
  """
  Identifies and prints the indices of dragonfly species expected to have
  reduced pterostigmata based on their migratory and gliding ecology.
  """

  # Species known for long-distance migration and a gliding flight style
  # are expected to have reduced pterostigmata to improve aerodynamic efficiency.
  # The indices below correspond to:
  # 3: Macrodiplax balteata
  # 4: Pantala flavescens
  # 8: Sympetrum corruptum
  # 10: Tholymis tillarga
  
  migratory_glider_indices = [3, 4, 8, 10]
  
  # The prompt asks for the output as a comma-separated string of indices.
  # The following line converts each number in the list to a string and joins them with a comma.
  output_string = ",".join(map(str, migratory_glider_indices))
  
  print("The indices of species expected to have reduced pterostigmata are:")
  print(output_string)

get_species_indices()