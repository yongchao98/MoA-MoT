def solve_identities():
    """
    This function determines the correct identities based on the mathematical derivation.
    The derivation shows that identities 4, 6, 7, 8, 10, 11, and 12 are the ones that
    necessarily follow from the given premises.
    """
    # List of numbers for the identities that are necessarily true.
    true_identity_numbers = [4, 6, 7, 8, 10, 11, 12]

    # The problem asks for the answer as a comma-separated string of these numbers.
    # The numbers should be in increasing order and without spaces.
    answer_string = ",".join(map(str, sorted(true_identity_numbers)))

    # The prompt says: "remember in the final code you still need to output each number in the final equation!"
    # I interpret this as printing the final comma-separated string, which represents the solution set.
    print(answer_string)

solve_identities()