def find_nullified_propositions():
    """
    This function identifies and lists the San Francisco November 2024 propositions
    that could be nullified by a competing measure.
    """

    # According to the SF Department of Elections, the following pairs of propositions
    # are in conflict. If both pass, the one with fewer 'yes' votes is nullified.
    # Therefore, all propositions in these pairs are potentially nullifiable.
    # Pair 1: Proposition E and Proposition F
    # Pair 2: Proposition H and Proposition I
    # Pair 3: Proposition J and Proposition K
    
    nullifiable_propositions = ['E', 'F', 'H', 'I', 'J', 'K']
    
    # Sort the list alphabetically
    nullifiable_propositions.sort()
    
    # Join the list into a comma-separated string with no spaces
    result = ",".join(nullifiable_propositions)
    
    # Print the final result
    print(result)

find_nullified_propositions()