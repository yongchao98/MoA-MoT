def find_nullifiable_propositions():
    """
    Identifies and lists the San Francisco Nov 2024 propositions
    that could be nullified by competing measures.
    """
    # According to the SF Department of Elections Voter Guide,
    # the following pairs of propositions are competing measures.
    # If both pass, the one with fewer 'yes' votes is nullified.
    # Therefore, all propositions in these pairs could potentially be nullified.
    conflicting_pairs = [
        ('C', 'D'),  # Real Property Transfer Tax
        ('F', 'G'),  # Business Tax Reform
        ('J', 'K')   # Supervisorial District Map
    ]

    # Flatten the list of pairs into a single list of propositions
    nullifiable_props = []
    for pair in conflicting_pairs:
        nullifiable_props.extend(pair)

    # Sort the list alphabetically
    nullifiable_props.sort()

    # Join the list into a comma-separated string with no spaces
    result_string = ",".join(nullifiable_props)

    print(result_string)

find_nullifiable_propositions()