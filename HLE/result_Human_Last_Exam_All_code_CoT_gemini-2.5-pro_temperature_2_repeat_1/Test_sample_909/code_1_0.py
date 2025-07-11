def solve_and_print_electric_field():
  """
  This function prints the derived expressions for the electric field
  in each region of the resistor, matching the format of the correct answer choice.
  """
  
  # Expression for the electric field in Region 1 (0 < phi < pi/2)
  E1_numerator = "2 * sigma_2 * V_0"
  E1_denominator = "r * pi * (sigma_1 + sigma_2)"
  E1_direction = "i_phi"
  
  print("Electric Field in Region 1:")
  print(f"E_1 = ({E1_numerator}) / ({E1_denominator}) * {E1_direction}")
  print("-" * 30)

  # Expression for the electric field in Region 2 (pi/2 < phi < pi)
  E2_numerator = "2 * sigma_1 * V_0"
  E2_denominator = "r * pi * (sigma_1 + sigma_2)"
  E2_direction = "i_phi"

  print("Electric Field in Region 2:")
  print(f"E_2 = ({E2_numerator}) / ({E2_denominator}) * {E2_direction}")

solve_and_print_electric_field()