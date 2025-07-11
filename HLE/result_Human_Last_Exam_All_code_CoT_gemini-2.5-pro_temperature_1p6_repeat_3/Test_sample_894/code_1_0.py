def find_nullified_propositions():
    """
    Identifies and lists the San Francisco November 2024 propositions
    that would be nullified if a competing proposition received more votes.
    """
    # Based on research of the SF November 2024 ballot, the following propositions
    # have conflict clauses with another proposition on the same ballot.
    # If both pass, the one with fewer 'yes' votes is nullified.
    # The pairs are (E, F), (G, H), and (M, N).
    # Therefore, any of these could be the one with fewer votes and be nullified.
    
    nullifiable_propositions = ['E', 'F', 'G', 'H', 'M', 'N']
    
    # Sort the list alphabetically to ensure correct order
    nullifiable_propositions.sort()
    
    # Join the list into a comma-separated string with no spaces
    result_string = ",".join(nullifiable_propositions)
    
    print(result_string)

find_nullified_propositions()