def calculate_fz(z):
  """
  Calculates the value of the pdf f_Z(z) for a given z.
  The formula is f_Z(z) = 1.5 - 3*z + 3*z**2.
  """
  term1 = 1.5
  term2 = 3 * z
  term3 = 3 * (z ** 2)
  
  result = term1 - term2 + term3
  
  print(f"The probability density function is f_Z(z) = 1.5 - 3*z + 3*z^2")
  print(f"For z = {z}:")
  print(f"f_Z({z}) = {term1} - {term2} + {term3}")
  print(f"f_Z({z}) = {result}")

# The point at which to evaluate the pdf
z_value = 0.2
calculate_fz(z_value)