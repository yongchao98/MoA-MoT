import itertools

def solve_ita_matrix_code():
    """
    This function determines the correct ITA Matrix extension code based on a set of rules.
    """
    # Step 1: Define the components based on the query.
    # - Direct flight (no connections on outbound): 'X'
    # - Operating carrier is Delta: 'O:DL'
    # - Eligible fare classes for Platinum Companion Certificate: 'L','U','T','X','V'
    #   The ITA syntax is '/f bc=FARE1|FARE2|...'
    
    # We apply lowercasing where possible as per the instructions.
    # 'X' becomes 'x', 'O:' becomes 'o:'.
    # Airline and fare class codes must remain uppercase.
    components = [
        'x',
        'o:DL',
        '/f bc=L|U|T|X|V'
    ]

    # Step 2: Generate all permutations of the components.
    # This ensures we can find the arrangement that meets the lexicographical requirement.
    all_permutations = list(itertools.permutations(components))

    # Step 3: Create the candidate strings from the permutations.
    # We join with a space, as it's the shortest separator. All resulting strings
    # will have the same minimal length.
    candidate_strings = [" ".join(p) for p in all_permutations]

    # Step 4: Sort the strings to find the lexicographically highest one.
    # Using key=str.lower ensures case-insensitive sorting as requested.
    candidate_strings.sort(key=str.lower)

    # The highest value is the last one in the sorted list.
    final_code = candidate_strings[-1]

    # Step 5: Print the final result clearly.
    print("The ITA Matrix outbound extension code for direct Delta flights eligible for Platinum Companion Certificates is constructed from the following components:")
    for component in components:
        print(f"- {component}")
    
    print("\nTo satisfy the conditions of minimal length and highest case-insensitive lexicographical value, we find the best ordering of these components.")
    print("\nThe optimal extension code is:")
    print(final_code)

if __name__ == '__main__':
    solve_ita_matrix_code()