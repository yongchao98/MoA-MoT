def find_nullifiable_propositions():
    """
    This function identifies and lists the San Francisco November 2024 propositions
    that would be nullified if a competing measure receives more votes.

    Based on the City Attorney's digests:
    - Proposition K and Proposition L are competing tax measures.
    - Proposition M and Proposition N are competing housing measures.

    In both cases, if both propositions in a pair pass, the one with more "yes" votes
    goes into effect, nullifying the other. Therefore, K, L, M, and N are all
    propositions that could be nullified.
    """
    
    # List of the letters for the propositions that could be nullified
    nullifiable_propositions = ['K', 'L', 'M', 'N']
    
    # Ensure the list is in alphabetical order
    nullifiable_propositions.sort()
    
    # Join the list into a single string, separated by commas with no spaces
    result_string = ",".join(nullifiable_propositions)
    
    # Print the final result
    print(result_string)

find_nullifiable_propositions()