def get_conflicting_propositions():
    """
    Identifies and prints the letters of the San Francisco November 2024 propositions
    that could be nullified if a competing measure passes with more votes.
    """
    # According to the SF Department of Elections, the following propositions have
    # clauses that nullify them if a competing measure receives more votes:
    # Prop C vs. Prop D
    # Prop I vs. Prop J
    # Prop L vs. Prop M
    # Therefore, all six are capable of being nullified.
    nullifiable_propositions = ['C', 'D', 'I', 'J', 'L', 'M']

    # Sort the list alphabetically to ensure consistent ordering.
    nullifiable_propositions.sort()

    # Join the list elements into a single string, separated by commas, with no spaces.
    result_string = ",".join(nullifiable_propositions)

    # Print the final result.
    print(result_string)

get_conflicting_propositions()