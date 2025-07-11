def get_jzdec_formula():
  """
  This function returns the Linear Logic formula for F(r, z, nz).
  """
  # S_z represents the state to jump to if the counter is zero.
  # C_r represents the counter to be tested.
  # S_nz represents the state to jump to if the counter is non-zero.
  # '&' is the additive conjunction 'with'.
  # '⊸' is the linear implication 'lollipop'.
  formula = "S_z & (C_r ⊸ S_nz)"
  
  # The problem asks to output each number in the final equation.
  # In this context, 'r', 'z', and 'nz' are the symbolic parameters.
  # We will print the formula with these symbols.
  print("The formula F(r, z, nz) is:")
  print(formula)

get_jzdec_formula()