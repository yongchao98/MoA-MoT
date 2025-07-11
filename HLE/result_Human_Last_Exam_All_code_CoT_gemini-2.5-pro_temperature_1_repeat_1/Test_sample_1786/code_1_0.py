import sys

# Define unicode characters for logical connectives for better readability
TENSOR = "\u2297"  # ⊗
WITH = "&"         # &
LOLLI = "\u22b8"   # ⊸

def generate_jzdec_formula(r, z, nz):
    """
    Constructs and prints the Linear Logic formula F(r, z, nz)
    for a JZDEC instruction, explaining its components.

    The Minsky instruction l: JZDEC(r, z, nz) is encoded as:
    !(S_l ⊸ F(r, z, nz))

    This function determines and prints the formula for F(r, z, nz).

    Args:
        r (int or str): The register index to test.
        z (int or str): The state label for the zero-counter case.
        nz (int or str): The state label for the non-zero-counter case.
    """

    # The formula uses the additive conjunction '&' (with) to create an
    # internal choice for the proof, representing the conditional nature
    # of the JZDEC instruction. The context (i.e., the value of counter r)
    # determines which of the two branches can lead to a successful proof.

    # 1. The Zero Branch: (S_z ⊗ Z_r)
    # This branch corresponds to the case where counter r is 0.
    # It produces the new state literal S_z and the zero-test literal Z_r.
    # The axioms in Delta ensure that Z_r can only be eliminated (allowing the
    # proof to proceed) if there are no C_r literals in the context.
    # If C_r is present, this branch of the proof fails.
    zero_branch = f"(S_{z} {TENSOR} Z_{r})"

    # 2. The Non-Zero Branch: (C_r ⊸ S_nz)
    # This branch corresponds to the case where counter r is positive.
    # The linear implication '⊸' consumes one C_r literal (decrementing the
    # counter) and produces the new state literal S_nz.
    # If there are no C_r literals, this branch cannot be used and fails.
    nonzero_branch = f"(C_{r} {LOLLI} S_{nz})"

    # The complete formula F(r, z, nz) combines these two branches.
    # A valid Minsky machine computation path will correspond to a proof
    # path where the correct branch of the '&' is chosen at each step.
    formula = f"{zero_branch} {WITH} {nonzero_branch}"

    print("The formula F(r, z, nz) that encodes the JZDEC instruction is:")
    print(formula)
    # The following line explicitly prints all the numbers used in the formula,
    # as requested by the prompt.
    print(f"\nThe numbers in this equation are: {z}, {r}, {r}, {nz}")


# --- Example Usage ---
# Suppose we want to encode the instruction: l: JZDEC(r=2, z=4, nz=5)
# Here, the register is 2.
# The 'zero' branch jumps to state 4.
# The 'non-zero' branch jumps to state 5.
if __name__ == '__main__':
    # You can change these values to see the formula for different instructions.
    register_to_test = 2
    zero_label = 4
    nonzero_label = 5

    generate_jzdec_formula(register_to_test, zero_label, nonzero_label)