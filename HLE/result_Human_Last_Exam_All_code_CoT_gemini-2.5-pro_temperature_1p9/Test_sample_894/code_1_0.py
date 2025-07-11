def find_nullified_propositions():
    """
    Identifies and prints the letters of the November 2024 San Francisco
    propositions that could be nullified by a competing measure.

    - Proposition E and F are competing measures; the one with fewer votes is nullified.
    - Proposition J is nullified if Proposition C passes with more votes.
    """
    # List of propositions that could be nullified
    nullified_propositions = ['E', 'F', 'J']

    # Sort them alphabetically
    nullified_propositions.sort()

    # Join them into a comma-separated string with no spaces
    result_string = ",".join(nullified_propositions)

    # Print the final result
    print(result_string)

find_nullified_propositions()