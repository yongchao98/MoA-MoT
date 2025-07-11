def find_nullifiable_propositions():
    """
    Identifies and lists San Francisco propositions that could be nullified by competing measures
    in the November 2024 election.
    """
    # This dictionary represents the pairs of competing propositions.
    # If Prop A competes with Prop B, both A and B can be nullified.
    # Data is based on the official SF Department of Elections Voter Information Pamphlet.
    competing_propositions = {
        "D": "N",
        "I": "J",
        "K": "M",
    }

    nullifiable = set()

    # If a proposition competes with another, both are added to the set.
    # Using a set automatically handles any duplicates.
    for prop1, prop2 in competing_propositions.items():
        nullifiable.add(prop1)
        nullifiable.add(prop2)

    # Convert the set to a list and sort it alphabetically.
    sorted_nullifiable_list = sorted(list(nullifiable))

    # Join the list into a comma-separated string with no spaces.
    result = ",".join(sorted_nullifiable_list)

    print(result)

find_nullifiable_propositions()