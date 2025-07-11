def solve():
    """
    Solves the problem by finding all maximal subsets of properties and ordering them.
    """

    # List of maximal subsets found through mathematical analysis.
    # Each subset is represented as a list of characters.
    maximal_subsets = [
        ['A', 'B'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['A', 'D'],
        ['A', 'E'],
        ['B', 'D'],
        ['C', 'D']
    ]

    # The properties within each subset are already ordered alphabetically.
    # Now, sort the list of subsets lexicographically.
    # Python's default list sorting is lexicographical.
    sorted_subsets = sorted(maximal_subsets)

    # Format the output string as per the example.
    # e.g., "{A,B}, {A,C,D}"
    formatted_parts = []
    for subset in sorted_subsets:
        # Join elements inside a set with ", " and enclose in braces
        formatted_set = "{" + ", ".join(subset) + "}"
        formatted_parts.append(formatted_set)
    
    # Join the formatted sets with ", "
    final_answer = ", ".join(formatted_parts)
    
    print(final_answer)

solve()