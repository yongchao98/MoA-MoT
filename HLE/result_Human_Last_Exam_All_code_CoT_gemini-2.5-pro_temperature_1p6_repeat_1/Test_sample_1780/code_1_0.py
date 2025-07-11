import itertools

def count_valid_arrangements():
    """
    Calculates the number of ways the letters "L", "N", "S", and "W" can be
    arranged based on a connection rule derived from their shapes.
    """
    # Step 1: Define letters and their tip properties.
    # The tuple represents (leftmost_tip_position, rightmost_tip_position).
    # 'T' stands for Top, 'B' for Bottom.
    letter_properties = {
        'L': ('T', 'B'),
        'N': ('T', 'B'),
        'S': ('T', 'B'),
        'W': ('T', 'T')
    }

    # The connection rule: The right tip of the first letter must have a
    # different position than the left tip of the second letter.
    def is_valid_connection(char1, char2):
        right_tip_1 = letter_properties[char1][1]
        left_tip_2 = letter_properties[char2][0]
        return right_tip_1 != left_tip_2

    # Step 2: Generate all permutations and check their validity.
    letters = ['L', 'N', 'S', 'W']
    all_permutations = itertools.permutations(letters)
    
    valid_arrangement_count = 0
    
    for arrangement in all_permutations:
        # An arrangement is a tuple like ('L', 'N', 'S', 'W')
        is_sequence_valid = True
        
        # Check connections between all adjacent letters in the arrangement.
        for i in range(len(arrangement) - 1):
            char1 = arrangement[i]
            char2 = arrangement[i+1]
            if not is_valid_connection(char1, char2):
                is_sequence_valid = False
                break
        
        # If all connections in the sequence are valid, increment the counter.
        if is_sequence_valid:
            valid_arrangement_count += 1
            
    # Step 3: Print the final count. The problem asks for the number of arrangements.
    # As per our logical deduction, this problem simplifies to finding the number
    # of permutations of {L, N, S}, which is 3!, as 'W' must always be last.
    # 3! = 3 * 2 * 1 = 6. The code calculates this by checking all possibilities.
    print(f"The number of ways is {valid_arrangement_count}.")

count_valid_arrangements()