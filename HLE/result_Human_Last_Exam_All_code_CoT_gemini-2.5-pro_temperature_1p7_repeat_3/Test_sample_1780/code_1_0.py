import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, W can be arranged
    based on a specific connection rule.
    """
    letters = ['L', 'N', 'S', 'W']

    # Adjacency list representing allowed connections (A -> B)
    # Based on the rule: S->W is valid, W->S is invalid.
    # This implies:
    # 1. Letters ending at the bottom (L, N, S) can connect to any letter.
    # 2. The letter ending at the top (W) can connect to any letter except S.
    adj = {
        'L': ['N', 'S', 'W'],
        'N': ['L', 'S', 'W'],
        'S': ['L', 'N', 'W'],
        'W': ['L', 'N']
    }

    # Generate all permutations of the four letters
    all_permutations = list(itertools.permutations(letters))
    
    valid_arrangements = []
    
    # Check each permutation for validity
    for p in all_permutations:
        p1, p2, p3, p4 = p
        # A permutation is valid if each consecutive pair is a valid connection
        if p2 in adj[p1] and p3 in adj[p2] and p4 in adj[p3]:
            valid_arrangements.append("".join(p))
            
    # Count how many valid arrangements start with each letter
    counts_by_start_letter = {'L': 0, 'N': 0, 'S': 0, 'W': 0}
    for arr in valid_arrangements:
        counts_by_start_letter[arr[0]] += 1
        
    total_count = len(valid_arrangements)
    
    # Construct the final equation string
    parts = []
    for letter in sorted(counts_by_start_letter.keys()):
        count = counts_by_start_letter[letter]
        if count > 0:
            parts.append(str(count))
            
    equation_str = " + ".join(parts) + f" = {total_count}"
    
    print(f"Number of arrangements starting with L: {counts_by_start_letter['L']}")
    print(f"Number of arrangements starting with N: {counts_by_start_letter['N']}")
    print(f"Number of arrangements starting with S: {counts_by_start_letter['S']}")
    print(f"Number of arrangements starting with W: {counts_by_start_letter['W']}")
    print("\nThe final equation is:")
    print(equation_str)
    
    print("\nThe total number of ways is:")
    print(total_count)


solve_letter_arrangement()