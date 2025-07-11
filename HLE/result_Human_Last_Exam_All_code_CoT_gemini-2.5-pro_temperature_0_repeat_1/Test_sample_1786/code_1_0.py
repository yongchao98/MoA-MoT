def generate_jzdec_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for a JZDEC instruction.

    Args:
        r (int): The register index (1-based).
        z (str): The label for the zero case.
        nz (str): The label for the non-zero case.
    """
    # The formula for the case where counter r is zero.
    # It produces the new state literal S_z.
    zero_case_formula = f"S_{z}"

    # The formula for the case where counter r is non-zero.
    # It consumes a counter resource C_r and produces the new state literal S_nz.
    # C_r ⊸ S_nz means "if you give me a C_r, I will give you an S_nz".
    nonzero_case_formula = f"(C_{r} ⊸ S_{nz})"

    # The full formula for F uses the '&' (with) connective to represent
    # the choice between the two cases. The context (presence or absence of C_r)
    # determines which branch can lead to a successful proof.
    final_formula = f"{zero_case_formula} & {nonzero_case_formula}"

    print("The formula F(r, z, nz) is:")
    print(final_formula)

if __name__ == '__main__':
    # Example usage of the function with concrete values.
    # Let's model the instruction l: JZDEC(r=1, z='halt', nz='loop')
    register_r = 1
    state_z = 'halt'
    state_nz = 'loop'
    
    # The problem asks to output each number in the final equation.
    # Here, the number is the register index 'r'.
    print(f"For the JZDEC instruction on register r = {register_r}:")
    print(f"If counter {register_r} is zero, go to state '{state_z}'.")
    print(f"If counter {register_r} is non-zero, go to state '{state_nz}'.")
    print("-" * 20)
    
    generate_jzdec_formula(register_r, state_z, state_nz)
