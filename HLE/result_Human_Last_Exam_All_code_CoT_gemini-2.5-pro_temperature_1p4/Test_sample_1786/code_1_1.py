import collections

def print_formula(r, z, nz):
    """
    Constructs and prints the Linear Logic formula F(r, z, nz).

    Args:
      r: The register index (integer).
      z: The state label for the zero case (string).
      nz: The state label for the non-zero case (string).
    """
    
    # Using unicode for logic symbols for better readability
    tensor = "\u2297"  # ⊗
    multimap = "\u22B8" # ⊸
    with_op = "&"      # &

    # The non-zero branch of the conditional.
    # It consumes a C_r resource (decrementing the counter) and produces the new state S_nz.
    non_zero_branch = f"C_{r} {multimap} S_{nz}"

    # The zero branch of the conditional.
    # It produces the new state S_z and a zero-test token Z_r.
    # This branch is only provable if no C_r resource is present in the context.
    zero_branch = f"S_{z} {tensor} Z_{r}"

    # The full formula F combines both branches with the additive conjunction '&'.
    # This creates an internal choice, forcing the selection of the branch
    # that is provable in the current context (i.e., based on whether C_r is present).
    final_formula = f"({non_zero_branch}) {with_op} ({zero_branch})"

    print(f"The formula F({r}, {z}, {nz}) is:")
    print(final_formula)
    print("\nIn this notation:")
    print(f"*   S_{z} represents the state for the zero case.")
    print(f"*   S_{nz} represents the state for the non-zero case.")
    print(f"*   C_{r} represents one unit of counter {r}.")
    print(f"*   Z_{r} is a special token to test if counter {r} is zero.")
    print(f"*   '{multimap}' is the linear implication (consumes resource on the left, produces on the right).")
    print(f"*   '{tensor}' is the multiplicative conjunction (puts both resources in the context).")
    print(f"*   '&' is the additive conjunction (provides a choice between the two formulas).")
    

# The specific components of the instruction from the prompt.
register_index = 'r'
zero_state_label = 'z'
non_zero_state_label = 'nz'

print_formula(register_index, zero_state_label, non_zero_state_label)
