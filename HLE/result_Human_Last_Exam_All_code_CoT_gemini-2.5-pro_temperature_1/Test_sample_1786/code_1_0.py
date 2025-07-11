def generate_formula(r, z, nz):
  """
  Generates the Linear Logic formula F(r, z, nz) for the JZDEC instruction.

  Args:
    r: The register index (integer).
    z: The state literal for the zero case (string).
    nz: The state literal for the non-zero case (string).
  """
  # The formula is a conjunction of two cases, represented by '&' (with).
  
  # Case 1: The counter x_r is zero.
  # We produce the new state S_z and a special atom Z_r for the zero test.
  # S_z stands for the new instruction register z.
  # Z_r is used to verify that counter r is indeed zero.
  zero_case = f"(S_{z} \u2297 Z_{r})" # Using unicode for tensor: ⊗
  
  # Case 2: The counter x_r is non-zero.
  # We must consume a C_r resource (decrement) and produce the new state S_nz.
  # This is modeled by a linear implication.
  # C_r stands for one unit in counter r.
  # S_nz stands for the new instruction register nz.
  nonzero_case = f"(C_{r} \u22b8 S_{nz})" # Using unicode for multimap: ⊸
  
  # The & connective provides an internal choice for the proof. The proof will only
  # succeed if it picks the branch that matches the state of the counter.
  # If x_r > 0, choosing the zero_case branch leaves C_r and Z_r as unremovable
  # resources, causing the proof to fail.
  # If x_r = 0, choosing the nonzero_case branch is impossible as the C_r resource
  # is not available.
  final_formula = f"{zero_case} & {nonzero_case}"
  
  print("The formula F(r, z, nz) is:")
  print(final_formula)

# Example usage with placeholder values r=1, z='z', nz='nz'
generate_formula(1, 'z', 'nz')