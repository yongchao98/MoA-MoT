def generate_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for the JZDEC instruction.

    Args:
        r (int or str): The register index (1-based).
        z (str): The state label for the zero case.
        nz (str): The state label for the non-zero case.
    """
    
    # Formula for the zero-branch:
    # If counter r is zero, produce the new state S_z and the zero-test literal Z_r.
    # The Z_r literal, using the axioms in Delta, will ensure this branch only succeeds
    # if no C_r literals are present, and will get stuck otherwise.
    F_zero = f"(S_{z} ⊗ Z_{r})"

    # Formula for the non-zero-branch:
    # This is a linear implication: given a C_r, produce an S_nz.
    # This means it consumes one C_r literal (decrement) and transitions to state nz.
    # This branch will get stuck if no C_r literal is available.
    F_nonzero = f"(C_{r} ⊸ S_{{nz}})"
    
    # The full formula uses the additive conjunction '&' (With) to combine the two
    # branches. This represents an external choice: the context of the proof
    # (i.e., whether C_r is present or not) determines which branch can succeed.
    F_final = f"{F_zero} & {F_nonzero}"
    
    print(f"The formula F({r},{z},{nz}) that encodes the JZDEC instruction is:")
    print(F_final)

# Example usage:
# Let's find the formula for an instruction l: JZDEC(1, "state_z", "state_nz")
r_val = 1
z_val = "z"
nz_val = "nz"
generate_formula(r_val, z_val, nz_val)

print("\nExplanation:")
print(f"F({r_val},{z_val},{nz_val}) = (S_{z_val} ⊗ Z_{r_val}) & (C_{r_val} ⊸ S_{{nz_val}})")
print("1. S_l is consumed and F is produced. The logic now offers a choice between the two branches of the '&'.")
print(f"2. Zero Case (x_{r_val} = 0):")
print(f"   - The right branch `(C_{r_val} ⊸ S_{{nz_val}})` gets stuck as no C_{r_val} resource is available.")
print(f"   - The left branch `(S_{z_val} ⊗ Z_{r_val})` succeeds. Z_{r_val} checks for absence of C_{r_val} and consumes other counter literals.")
print(f"3. Non-Zero Case (x_{r_val} > 0):")
print(f"   - The left branch `(S_{z_val} ⊗ Z_{r_val})` gets stuck because Z_{r_val} and C_{r_val} cannot be reduced together.")
print(f"   - The right branch `(C_{r_val} ⊸ S_{{nz_val}})` succeeds, consuming one C_{r_val} (decrement) and producing S_{{nz_val}}.")

<<<The formula F(r,z,nz) that encodes the JZDEC instruction is:
(S_z ⊗ Z_r) & (C_r ⊸ S_{nz})>>>