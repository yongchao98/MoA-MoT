import itertools

def solve_letter_arrangement():
    """
    Finds the number of valid arrangements for the letters L, N, S, W based on tip connection rules.
    """

    # Step 1 & 2: Define letter properties based on vertical alignment of their tips.
    # 'left' refers to the leftmost tip of the letter, 'right' to the rightmost tip.
    # The value represents the vertical position of that tip.
    letter_tips = {
        'L': {'left': 'Bottom', 'right': 'Middle'},
        'N': {'left': 'Top', 'right': 'Bottom'},
        'S': {'left': 'Bottom', 'right': 'Top'},
        'W': {'left': 'Top', 'right': 'Top'}
    }

    letters = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Step 3: Generate all possible arrangements (permutations) of the letters.
    all_arrangements = list(itertools.permutations(letters))

    # Step 4: Check each arrangement for valid connections.
    for arr in all_arrangements:
        # An arrangement is a sequence like (C1, C2, C3, C4)
        c1, c2, c3, c4 = arr[0], arr[1], arr[2], arr[3]

        # Check if C1 can connect to C2
        conn1_to_2 = letter_tips[c1]['right'] == letter_tips[c2]['left']
        # Check if C2 can connect to C3
        conn2_to_3 = letter_tips[c2]['right'] == letter_tips[c3]['left']
        # Check if C3 can connect to C4
        conn3_to_4 = letter_tips[c3]['right'] == letter_tips[c4]['left']

        # If all connections in the chain are valid, increment the counter.
        if conn1_to_2 and conn2_to_3 and conn3_to_4:
            valid_arrangement_count += 1

    # Step 5: Print the final count.
    print(valid_arrangement_count)

solve_letter_arrangement()