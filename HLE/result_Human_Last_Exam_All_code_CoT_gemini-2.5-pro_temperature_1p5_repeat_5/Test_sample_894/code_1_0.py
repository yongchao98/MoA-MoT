def find_nullifiable_propositions():
    """
    Identifies and lists San Francisco propositions that could be nullified by competing measures
    in the November 2024 election.
    """
    # According to the official SF Voter Information Pamphlet, Propositions B and E are competing measures.
    # If both pass, the one with more 'yes' votes becomes law, and the other is nullified.
    # This means both 'B' and 'E' are potentially nullifiable.
    competing_pairs = [
        ('B', 'E')
    ]

    # Use a set to store the unique propositions that can be nullified.
    nullifiable_props = set()

    # Add both propositions from each competing pair to the set.
    for prop1, prop2 in competing_pairs:
        nullifiable_props.add(prop1)
        nullifiable_props.add(prop2)

    # Sort the propositions alphabetically for the final output.
    sorted_props = sorted(list(nullifiable_props))

    # Format the list as a comma-separated string with no spaces.
    result_string = ",".join(sorted_props)

    print("List of propositions that could be nullified, in alphabetical order:")
    print(result_string)

if __name__ == "__main__":
    find_nullifiable_propositions()