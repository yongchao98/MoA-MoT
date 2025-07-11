def find_nullified_propositions():
    """
    Identifies and prints the San Francisco November 2024 propositions
    that would be nullified if a competing measure received more votes.
    """
    # Based on research, the following propositions are competing measures.
    # If both propositions in a pair pass, the one with fewer 'yes'
    # votes is nullified. Therefore, each of these could potentially be nullified.
    # Competing Pair 1: Prop H vs Prop I
    # Competing Pair 2: Prop K vs Prop L
    potentially_nullified_letters = ['H', 'I', 'K', 'L']

    # Sort the letters alphabetically, although they are already in order.
    potentially_nullified_letters.sort()

    # Join the list into a comma-separated string with no spaces.
    output_string = ",".join(potentially_nullified_letters)

    print(output_string)

find_nullified_propositions()