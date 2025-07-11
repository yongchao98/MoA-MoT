def identify_reagents():
  """
  This function identifies and prints the reagents A and B from the provided
  chemical reaction scheme.
  """

  # Reagent A introduces the N-NH2 group.
  reagent_A_name = "Hydrazine"
  reagent_A_formula = "H2N-NH2"

  # Reagent B introduces the N-propyl group and facilitates ring formation.
  reagent_B_name = "n-Propylamine"
  reagent_B_formula = "CH3CH2CH2NH2"

  print(f"Reagent A is: {reagent_A_name} (Formula: {reagent_A_formula})")
  print(f"Reagent B is: {reagent_B_name} (Formula: {reagent_B_formula})")

identify_reagents()