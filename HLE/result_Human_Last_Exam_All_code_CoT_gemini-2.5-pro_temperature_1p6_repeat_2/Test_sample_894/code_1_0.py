def find_nullified_propositions():
    """
    Identifies and prints the San Francisco Nov 2024 propositions
    that could be nullified by a competing measure.
    """
    # Based on the SF Department of Elections Voter Guide, the following pairs
    # of propositions are competing. If both propositions in a pair pass, the
    # one with fewer 'yes' votes is nullified.
    # Competing Pair 1: Proposition C vs. Proposition D
    # Competing Pair 2: Proposition H vs. Proposition I
    # Competing Pair 3: Proposition K vs. Proposition L
    # Therefore, any of these six propositions could be nullified.
    
    nullifiable_props = ['C', 'D', 'H', 'I', 'K', 'L']
    
    # Sort them alphabetically to ensure consistent ordering.
    nullifiable_props.sort()
    
    # Join the list into a comma-separated string with no spaces.
    result = ",".join(nullifiable_props)
    
    print(result)

find_nullified_propositions()