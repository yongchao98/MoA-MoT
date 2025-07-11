def solve_electron_energy():
  """
  This function calculates and prints the minimum energy of electron 1
  for the described process. The result is expressed in terms of the
  bandgap energy, E_g.
  """
  # The problem reduces to finding the minimum of E_1_initial, which is given by:
  # E_1_initial_min = 2*E_g + min(T_sum)
  # where min(T_sum) is the minimum possible sum of the kinetic energies of the
  # other particles involved, which is found to be E_g / 2.

  # Coefficients from the derivation
  coefficient_base = 2
  coefficient_additional = 0.5
  
  # Calculate the final coefficient
  total_coefficient = coefficient_base + coefficient_additional

  # Print the equation with all numbers
  print(f"The minimum initial energy for electron 1 (E_1_min) is derived from the conservation laws.")
  print(f"The calculation is as follows:")
  print(f"E_1_min = {coefficient_base} * E_g + {coefficient_additional} * E_g = {total_coefficient} * E_g")

solve_electron_energy()