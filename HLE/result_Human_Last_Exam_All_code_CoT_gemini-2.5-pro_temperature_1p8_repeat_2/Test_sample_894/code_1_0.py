import collections

def find_nullifiable_propositions():
    """
    Identifies and lists San Francisco propositions that could be nullified by a competing measure on the same ballot.

    The November 2024 San Francisco ballot includes propositions where passage is
    contingent on receiving more votes than a competing measure. This script
    identifies such propositions based on documented conflicts.
    """
    # According to the SF Department of Elections Voter Information Pamphlet,
    # Propositions H and J are competing measures. If both pass, the one with
    # more 'yes' votes goes into effect. Therefore, either one could be nullified.
    # We represent this relationship here. A value indicates a competitor.
    propositions_with_competitors = {
        'A': None,
        'B': None,
        'C': None,
        'D': None,
        'E': None,
        'F': None,
        'G': None,
        'H': 'J', # Competes with Proposition J
        'I': None,
        'J': 'H', # Competes with Proposition H
        'K': None,
        'L': None,
        'M': None,
        'N': None,
        'O': None
    }

    nullifiable_props = []
    # Iterate through the propositions and find any that have a competitor.
    for prop, competitor in propositions_with_competitors.items():
        if competitor is not None:
            nullifiable_props.append(prop)

    # Sort the list of nullifiable propositions alphabetically.
    nullifiable_props.sort()

    # Join the sorted list into a comma-separated string with no spaces.
    result = ",".join(nullifiable_props)
    
    print(f"The propositions that could be nullified are: {result}")

find_nullifiable_propositions()