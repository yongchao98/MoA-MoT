def is_filled_nilpotent(sylow_factors):
    """
    Checks if a finite nilpotent group is filled based on its Sylow factors.

    A finite nilpotent group G is the direct product of its Sylow p-subgroups.
    We represent G by a list of tuples (p, exp), where each tuple represents a
    Sylow p-subgroup P_p with its prime p and its exponent exp.

    Args:
        sylow_factors: A list of tuples, e.g., [(p1, exp1), (p2, exp2), ...].

    Returns:
        True if the group is filled, False otherwise.
    """
    num_factors = len(sylow_factors)

    if num_factors == 0:
        # The trivial group is filled.
        return True

    # Case 1: Check if G is a p-group of exponent p.
    # This means there is only one Sylow factor, and its exponent equals its prime.
    if num_factors == 1:
        p, exp = sylow_factors[0]
        if p == exp:
            return False  # This group is not filled.

    # Case 2: Check if G = P x Q, where P is a p-group of exponent p
    # and Q is a q-group of exponent q.
    # This means there are exactly two Sylow factors, and each has an exponent
    # equal to its prime.
    if num_factors == 2:
        p1, exp1 = sylow_factors[0]
        p2, exp2 = sylow_factors[1]
        if p1 == exp1 and p2 == exp2:
            return False # This group is not filled.

    # If the group does not match either of the non-filled cases, it is filled.
    return True

def main():
    """Main function to print the theorem and test example groups."""
    print("--- The Classification of Finite Filled Nilpotent Groups ---")
    print("""
A finite nilpotent group G is called 'filled' if every non-identity element
belongs to at least one maximal product-free set.

According to a theorem by Cairns and de Almeida (2021), a finite nilpotent
group G is filled if and only if it is NOT one of the following two types:

1. A p-group of exponent p (for a prime p).
   This is a group of prime-power order p^k where every non-identity
   element has order p.

2. A direct product G = P × Q, where P is a p-group of exponent p
   and Q is a q-group of exponent q (for distinct primes p and q).

Below are some examples tested using this classification.
A group is represented by its Sylow factors [(p_1, exp_1), (p_2, exp_2), ...].
""")

    # --- Dictionary of test groups ---
    # Key: Group name (for display)
    # Value: Sylow factor representation [(p, exp), ...]
    test_groups = {
        "C_5 (Cyclic group of order 5)": [(5, 5)],
        "C_2 x C_2 (Klein four-group)": [(2, 2)],
        "Q_8 (Quaternion group)": [(2, 4)],
        "C_6 (Cyclic group C_2 x C_3)": [(2, 2), (3, 3)],
        "C_12 (Cyclic group C_4 x C_3)": [(2, 4), (3, 3)],
        "C_30 (Cyclic group C_2 x C_3 x C_5)": [(2, 2), (3, 3), (5, 5)],
        "C_7 x Q_8": [(7, 7), (2, 4)]
    }

    # --- Test each group and print the result ---
    for name, factors in test_groups.items():
        is_filled = is_filled_nilpotent(factors)
        description = "Filled" if is_filled else "Not Filled"
        # Print the final equation for each group tested
        sylow_str = ' × '.join([f"P_{p}(exp={exp})" for p, exp in factors])
        print(f"Group: {name} ({sylow_str})")
        print(f"Result: This group is {description}.\n")

if __name__ == "__main__":
    main()
