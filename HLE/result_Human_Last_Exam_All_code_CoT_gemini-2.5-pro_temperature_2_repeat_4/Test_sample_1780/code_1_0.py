import itertools

def solve_letter_arrangement():
    """
    Calculates the number of valid arrangements of the letters L, N, S, W
    based on a specific tip-connection rule.
    """
    
    # Step 1: Define the properties of each letter's unconnected tips.
    # The tuple represents (left_tip_position, right_tip_position).
    # 't' stands for top, 'b' for bottom.
    tip_positions = {
        'L': ('t', 'b'),
        'N': ('b', 't'),
        'S': ('t', 'b'),
        'W': ('t', 't')
    }

    # Step 2 is implemented in this helper function.
    # The rule is that the right tip of the first letter and the left tip
    # of the second letter must have different vertical positions.
    def is_valid_connection(letter1, letter2):
        right_tip1 = tip_positions[letter1][1]
        left_tip2 = tip_positions[letter2][0]
        return right_tip1 != left_tip2

    # Step 3: Generate all possible arrangements (permutations) of the letters.
    letters = ["L", "N", "S", "W"]
    all_arrangements = itertools.permutations(letters)

    # Step 4: Count the arrangements that are valid.
    valid_arrangement_count = 0
    
    for arrangement in all_arrangements:
        # An arrangement (p1, p2, p3, p4) is valid if all three
        # consecutive connections (p1->p2, p2->p3, p3->p4) are valid.
        if (is_valid_connection(arrangement[0], arrangement[1]) and
            is_valid_connection(arrangement[1], arrangement[2]) and
            is_valid_connection(arrangement[2], arrangement[3])):
            valid_arrangement_count += 1
            
    # Step 5: Output the final number of valid arrangements.
    print(valid_arrangement_count)

solve_letter_arrangement()