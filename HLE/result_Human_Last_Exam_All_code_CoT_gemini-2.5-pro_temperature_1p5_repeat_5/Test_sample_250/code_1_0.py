def solve_scheme_properties():
    """
    Finds and prints the lexicographically sorted list of maximal subsets of properties
    for which a scheme exists.
    """
    # The maximal subsets found through mathematical reasoning.
    # Each subset is represented as a list of characters, sorted internally.
    maximal_subsets = [
        ['A', 'B'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['B', 'D']
    ]

    # Python's default list sort provides lexicographical ordering.
    maximal_subsets.sort()

    # Format the output string as a list of sets.
    # Example: "{A, B}, {A, C, D}"
    output_parts = []
    for subset in maximal_subsets:
        formatted_subset = "{" + ", ".join(subset) + "}"
        output_parts.append(formatted_subset)

    final_output = ", ".join(output_parts)
    print(final_output)

solve_scheme_properties()