def get_jzdec_formula(r, z, nz):
  """
  Constructs the Linear Logic formula F(r, z, nz) for a JZDEC instruction.

  Args:
    r: The index of the counter register (a variable).
    z: The label for the zero branch state (a variable).
    nz: The label for the non-zero branch state (a variable).
  
  Returns:
    A string representing the formula.
  """
  # S_z: Literal for the zero-branch state 'z'.
  # C_r: Literal for one unit in counter 'r'.
  # S_nz: Literal for the non-zero-branch state 'nz'.
  # &: The "with" connective (additive conjunction), providing a choice.
  # ⊸: The "lollipop" connective (linear implication), for resource transformation.
  return f"S_{z} & (C_{r} ⊸ S_{nz})"

# Define symbolic variables for the instruction parameters
register_variable = 'r'
zero_state_variable = 'z'
nonzero_state_variable = 'nz'

# Get the formula
formula = get_jzdec_formula(register_variable, zero_state_variable, nonzero_state_variable)

# The question asks for the formula F(r, z, nz) used in the axiom
# !(S_l ⊸ F(r, z, nz))
print("The appropriate formula F(r,z,nz) is:")
print(formula)

print("\nWhere:")
print(f"  r = {register_variable} (the index of the counter to test)")
print(f"  z = {zero_state_variable} (the label of the state to enter if counter r is zero)")
print(f"  nz = {nonzero_state_variable} (the label of the state to enter if counter r is non-zero)")

final_answer = 'S_z & (C_r ⊸ S_nz)'
print(f"\nFinal answer in required format: <<<{final_answer}>>>")