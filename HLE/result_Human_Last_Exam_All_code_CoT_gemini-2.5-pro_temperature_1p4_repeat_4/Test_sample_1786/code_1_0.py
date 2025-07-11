def generate_jzdec_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for a JZDEC instruction.

    Args:
        r (int): The register index (1-based).
        z (str): The state label for the zero case.
        nz (str): The state label for the non-zero case.

    Returns:
        str: The formatted Linear Logic formula.
    """
    # The formula structure is S_z & (C_r ⊸ S_nz)
    # Using Unicode for the linear implication arrow for clarity.
    formula = f"S_{z} & (C_{r} ⊸ S_{nz})"
    return formula

# --- Example Usage ---
# Let's model a specific JZDEC instruction from a hypothetical Minsky machine program.
# Instruction: l1: JZDEC(r=2, z='final', nz='l1_loop')

# Parameters for the instruction
register_index = 2
zero_case_label = 'final'
non_zero_case_label = 'l1_loop'

# Generate the corresponding formula F
final_formula = generate_jzdec_formula(register_index, zero_case_label, non_zero_case_label)

# Print the result, showing how the numbers and labels from the instruction
# are placed into the final formula.
print(f"For the Minsky machine instruction JZDEC(r={register_index}, z='{zero_case_label}', nz='{non_zero_case_label}'),")
print(f"the corresponding formula F({register_index}, '{zero_case_label}', '{non_zero_case_label}') is:")
print(final_formula)
