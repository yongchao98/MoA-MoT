def print_formula(r, z, nz):
  """
  Prints the Linear Logic formula for the JZDEC instruction.

  Args:
    r: The register index (integer).
    z: The label for the zero case (string).
    nz: The label for the non-zero case (string).
  """
  # The formula F(r, z, nz) needs to represent the two cases of the JZDEC instruction.
  # We use the additive conjunction '&' (with) to present a choice.
  
  # Case 1: The counter x_r is 0. The machine should jump to state z.
  # This is represented by the formula S_z, which provides the resource for the new state.
  zero_case = f"S_{z}"
  
  # Case 2: The counter x_r is greater than 0. The machine should decrement the counter
  # and jump to state nz. This operation of consuming a resource C_r and producing S_nz
  # is represented by the linear implication C_r --o S_nz.
  nonzero_case = f"(C_{r} --o S_{nz})"
  
  # The final formula combines these two choices using '&'.
  # When x_r = 0, there is no C_r resource, so the second choice is impossible,
  # forcing the first choice.
  # When x_r > 0, both choices are logically possible, but only the second one
  # corresponds to a valid Minsky machine step. The proof of undecidability
  # relies on the fact that choosing the "wrong" path will not lead to a valid halt.
  final_formula = f"{zero_case} & {nonzero_case}"
  
  print("The formula F(r,z,nz) that encodes the JZDEC instruction is:")
  print(final_formula)

# Example usage with placeholder values r=1, z='z', nz='nz'
print_formula(r='r', z='z', nz='nz')