import itertools

def solve_letter_arrangement():
    """
    Finds the number of valid arrangements for the letters L, N, S, W
    based on a specific tip-connection rule.
    """
    # Step 1 & 2: Define letter properties and the connection rule.
    # We use 0 for a 'low' tip and 1 for a 'high' tip.
    # The dictionary maps each letter to its left and right tip levels.
    letters = {
        'L': {'left': 0, 'right': 1}, # low -> high
        'N': {'left': 0, 'right': 1}, # low -> high
        'S': {'left': 1, 'right': 0}, # high -> low
        'W': {'left': 1, 'right': 1}  # high -> high
    }

    capitals = list(letters.keys())
    
    # Use a dictionary to store counts for arrangements starting with each letter.
    # The keys are sorted to ensure a deterministic output order.
    valid_counts_by_start = {key: 0 for key in sorted(capitals)}

    # Step 3: Generate and check all permutations.
    all_arrangements = itertools.permutations(capitals)

    for arr in all_arrangements:
        is_valid = True
        # An arrangement is a sequence (p1, p2, p3, p4).
        # We must check the connections p1->p2, p2->p3, and p3->p4.
        for i in range(len(arr) - 1):
            letter1 = arr[i]
            letter2 = arr[i+1]
            
            # Connection Rule: The right tip of letter1 must be at a different
            # level than the left tip of letter2.
            if letters[letter1]['right'] == letters[letter2]['left']:
                is_valid = False
                break
        
        if is_valid:
            # If the entire arrangement is valid, increment the counter
            # for its starting letter.
            start_letter = arr[0]
            valid_counts_by_start[start_letter] += 1
            
    # Step 4: Output the results in an "equation" format.
    # This shows the contribution from each possible starting letter.
    counts = list(valid_counts_by_start.values())
    total = sum(counts)
    
    # Creates the string "c1 + c2 + c3 + ..."
    equation_parts = " + ".join(map(str, counts))
    
    # We need to output each number in the final equation.
    # The order of numbers corresponds to the sorted letters: L, N, S, W.
    # valid_counts_by_start['L'] = 0
    # valid_counts_by_start['N'] = 0
    # valid_counts_by_start['S'] = 2
    # valid_counts_by_start['W'] = 0
    # The numbers in the equation are 0, 0, 2, and 0. The total is 2.
    print(f"Number of arrangements starting with L: {valid_counts_by_start['L']}")
    print(f"Number of arrangements starting with N: {valid_counts_by_start['N']}")
    print(f"Number of arrangements starting with S: {valid_counts_by_start['S']}")
    print(f"Number of arrangements starting with W: {valid_counts_by_start['W']}")
    print("---")
    print(f"Final Calculation: {equation_parts} = {total}")
    print(f"Total number of ways: {total}")


solve_letter_arrangement()