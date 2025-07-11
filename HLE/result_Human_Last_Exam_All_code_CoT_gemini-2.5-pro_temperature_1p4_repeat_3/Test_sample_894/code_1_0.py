def find_nullifiable_propositions():
    """
    Identifies and lists San Francisco propositions from the November 2024 election
    that could be nullified by other propositions.
    """
    # Define the groups of competing propositions. If propositions in a group pass,
    # the one with the most votes is enacted, nullifying the others in that group.
    competing_groups = [
        {'C', 'I', 'J'},  # Competing Real Estate Transfer Tax measures
        {'E', 'F'},      # Competing police policy and public safety measures
        {'K', 'L'}       # Competing housing development measures
    ]

    # Use a set to automatically handle duplicates as we collect all propositions
    nullifiable_props = set()

    for group in competing_groups:
        nullifiable_props.update(group)

    # Sort the proposition letters alphabetically
    sorted_props = sorted(list(nullifiable_props))

    # Join the sorted letters into a comma-separated string with no spaces
    result = ",".join(sorted_props)

    print(result)

find_nullifiable_propositions()