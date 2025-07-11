def get_splitting_pattern(n_neighbors):
  """
  Determines the splitting pattern based on the n+1 rule.
  """
  multiplicity = n_neighbors + 1
  if multiplicity == 1:
    return "singlet"
  elif multiplicity == 2:
    return "doublet"
  elif multiplicity == 3:
    return "triplet"
  elif multiplicity == 4:
    return "quartet"
  else:
    return f"{multiplicity}-peak multiplet"

def main():
  """
  Analyzes the structure of Compound 1 to determine the properties of the most deshielded proton.
  """
  # Step 1 & 2: Identify the most deshielded proton.
  # The most deshielded proton is the single H on the central aromatic ring,
  # located between the two nitrogen atoms. This is due to the strong deshielding
  # effect of the adjacent electronegative N atoms and the overall positive charge.
  most_deshielded_proton = "Proton on the central ring (meso-position)"

  # Step 3: Determine the integration.
  # There is only one such proton in the molecule.
  integration = 1

  # Step 4: Determine the splitting pattern.
  # The carbon atom this proton is attached to has no protons on its adjacent carbon atoms.
  # Therefore, the number of neighboring protons (n) is 0.
  num_neighboring_protons = 0
  splitting_pattern = get_splitting_pattern(num_neighboring_protons)

  # Print the final answer
  print(f"The highest deshielded proton peak in Compound 1 is a {splitting_pattern} with an integration of {integration}H.")

if __name__ == "__main__":
  main()