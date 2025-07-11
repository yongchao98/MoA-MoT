def find_nullified_propositions():
    """
    Identifies and prints the San Francisco 2024 propositions
    that could be nullified by a competing proposition.
    """
    # Propositions with conflict clauses where the one with more votes wins.
    # Source: San Francisco Department of Elections Voter Information Pamphlet.
    # Prop F conflicts with G.
    # Prop G conflicts with F.
    # Prop K conflicts with L.
    # Prop L conflicts with K.
    # Prop M conflicts with N.
    # Prop N conflicts with M.
    propositions = ['F', 'G', 'K', 'L', 'M', 'N']

    # The list is already in alphabetical order, but sorting ensures it.
    propositions.sort()

    # Join the list into a comma-separated string with no spaces.
    result = ",".join(propositions)

    print(result)

find_nullified_propositions()