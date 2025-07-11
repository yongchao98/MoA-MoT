import math

def get_density_constants():
  """
  This function describes the constants in the normalised invariant density function
  rho(x) = 1 / (2 * (1 - ln(2)) * (sqrt(x) + 1))
  for the map T(x) = 1/sqrt(x) mod 1.
  """
  
  # The final equation is rho(x) = 1 / (2 * (1 - ln(2)) * (sqrt(x) + 1))
  # The numbers in this equation are 1, 2, 1, 2, 1.
  
  print("The formula for the normalised density rho(x) is of the form:")
  print("1 / (c1 * (c2 - ln(c3)) * (sqrt(x) + c4))")
  
  print("\nThe values of the constants are:")
  print(1) # a
  print(2) # b
  print(1) # d
  print(2) # f
  print(1) # g

get_density_constants()
