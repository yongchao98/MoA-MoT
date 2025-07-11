import itertools

def solve():
    """
    Calculates the number of ways the letters L, N, S, W can be arranged
    based on specific connection rules.
    """
    
    # The connections are represented as a dictionary (an adjacency list for a graph).
    # The key is the starting letter, and the value is a list of letters it can connect to.
    # This graph is derived from analyzing the letter shapes and applying the specific
    # rule from the prompt: "S can connect to W but W cannot connect to S".
    connections = {
        'L': ['N'],
        'N': ['L', 'S', 'W'],
        'S': ['L', 'N', 'W'],
        'W': ['L', 'N']
    }

    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Generate all 4! = 24 possible permutations of the letters.
    all_permutations = list(itertools.permutations(letters))

    # Iterate through each permutation to see if it's a valid path.
    for p in all_permutations:
        # A permutation (p[0], p[1], p[2], p[3]) is a valid arrangement if
        # p[0]->p[1], p[1]->p[2], and p[2]->p[3] are all valid connections.
        try:
            if (p[1] in connections[p[0]] and
                p[2] in connections[p[1]] and
                p[3] in connections[p[2]]):
                valid_arrangement_count += 1
        except KeyError:
            # This handles cases where a letter might not have any outgoing connections,
            # though not applicable in this specific graph.
            continue
            
    # Print the final count of valid arrangements.
    print(valid_arrangement_count)

solve()