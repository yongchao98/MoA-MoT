import textwrap

def generate_jzdec_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for the JZDEC instruction.
    
    This function demonstrates the encoding of a Minsky machine's conditional
    branch instruction into a single Linear Logic formula, which is a key
    step in proving the undecidability of Linear Logic.

    Args:
        r (int): The index of the counter register to test (e.g., 1, 2, ...).
        z (str): The label of the instruction to jump to if the counter is zero.
        nz (str): The label of the instruction to jump to if the counter is non-zero.

    Returns:
        str: A string representing the complete formula F(r, z, nz).
    """

    # 1. The Non-Zero Case: (C_r -o S_nz)
    # This formula represents the action of decrementing counter r and moving to state nz.
    # The linear implication '-o' (multimap) means that to get S_nz,
    # we must consume a C_r resource. This is only possible if the counter r > 0.
    nonzero_case = f"(C_{r} \multimap S_{nz})"

    # 2. The Zero Case: (S_z tensor Z_r)
    # This formula represents the action of moving to state z after verifying counter r is zero.
    # It produces the new state S_z and a special "zero-test" literal Z_r.
    # The axioms provided in the problem for Z_r ensure it can only be resolved
    # to 1 (truth) if counter r is indeed zero. If r > 0, the proof branch gets stuck.
    # The symbol for the multiplicative conjunction 'tensor' is \otimes.
    zero_case = f"(S_{z} \otimes Z_{r})"

    # 3. Combine with Internal Choice '&' (with)
    # The final formula F combines these two mutually exclusive branches.
    # The '&' connective gives the prover a choice. However, the state of the
    # counters makes one of the choices lead to a dead end, forcing the
    # prover down the correct computational path.
    formula = f"{nonzero_case} & {zero_case}"
    return formula

def main():
    # Let's define a sample JZDEC instruction.
    # Instruction label 'l_test'
    # Counter register r = 2
    # Jump to 'l_final' if r is zero
    # Jump to 'l_loop' if r is non-zero
    r_val = 2
    z_label = "l_final"
    nz_label = "l_loop"

    print(textwrap.dedent(f"""
    Here is the construction of the formula F(r, z, nz) for the JZDEC instruction.
    
    Inputs:
      - Counter Register (r): {r_val}
      - Zero-branch State (z): {z_label}
      - Non-zero-branch State (nz): {nz_label}
    """))

    final_formula = generate_jzdec_formula(r_val, z_label, nz_label)

    print("The final formula F(r, z, nz) is:")
    print(final_formula)

if __name__ == "__main__":
    main()