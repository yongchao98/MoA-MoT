def get_nullifiable_propositions():
    """
    This function identifies and prints the letters of San Francisco propositions
    for the November 2024 election that could be nullified by a competing measure.

    According to the San Francisco Department of Elections, the following pairs
    of propositions are competing measures: (C, D) and (L, M). If both
    propositions in a pair pass, the one with fewer 'yes' votes is nullified.
    Therefore, all four are potentially nullifiable.
    """

    # List of propositions that can be nullified
    nullifiable_propositions = ['C', 'D', 'L', 'M']

    # Sort the list alphabetically
    nullifiable_propositions.sort()

    # Join the list into a comma-separated string with no spaces
    result_string = ",".join(nullifiable_propositions)

    # Print the final list
    print(result_string)

get_nullifiable_propositions()