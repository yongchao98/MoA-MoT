def generate_formula():
  """
  This function prints the formula for F(r, z, nz).
  """
  
  # Define the components of the formula symbolically
  # These represent the two cases of the JZDEC instruction.
  
  # Case 1: x_r > 0 (non-zero). Decrement counter r and go to state nz.
  # This is represented by consuming a C_r resource and producing an S_nz resource.
  nonzero_case = "C_r \multimap S_{nz}"
  
  # Case 2: x_r = 0 (zero). Go to state z.
  # This is represented by producing an S_z resource and initiating a zero-test
  # with the Z_r resource.
  zero_case = "S_z \otimes Z_r"
  
  # The full formula combines these two mutually exclusive cases using the
  # external choice connective '&' (with).
  formula = f"({nonzero_case}) & ({zero_case})"
  
  print("The formula for F(r, z, nz) is:")
  print(formula)

generate_formula()