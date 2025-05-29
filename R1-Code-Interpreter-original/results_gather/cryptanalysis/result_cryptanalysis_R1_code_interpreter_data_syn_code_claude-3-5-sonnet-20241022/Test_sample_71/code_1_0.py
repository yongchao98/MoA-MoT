def validate_against_all_conditions(candidate):
    conditions = [
        # Format: (guess, correct_nums, num_positions, correct_letters, letter_positions, num_size_info, letter_position_info)
        ("65HO", 1, "wrong", 0, [], "one_small", None),
        ("47IM", 0, [], 0, [], None, None),
        ("10GI", 0, [], 0, [], "both_small", "too_early"),
        ("70XG", 0, [], 0, [], None, None),
        ("12JU", 0, [], 1, [3], "both_small", "one_late"),
        ("52SZ", 0, [], 0, [], "both_small", "too_late"),
        ("64AV", 1, "wrong", 0, [], "one_small", None),
        ("85BO", 1, "correct", 0, [], "one_small", None),
        ("98RX", 1, "wrong", 0, [], "one_large", "too_late"),
        ("60JU", 1, "wrong", 1, [2], "one_small", "one_late"),
        ("75UC", 0, [], 0, [], None, None),
        ("79NC", 0, [], 1, "wrong", None, "one_early")
    ]
    
    for guess, correct_nums, num_pos, correct_letters, letter_pos, num_size, letter_info in conditions:
        # Check number of correct numbers
        correct_count = sum(1 for i in range(2) if guess[i] == candidate[i])
        wrong_pos_count = sum(1 for i in range(2) if guess[i] == candidate[1-i])
        
        if correct_nums == 1 and num_pos == "wrong" and wrong_pos_count != 1:
            return False
        if correct_nums == 1 and num_pos == "correct" and correct_count != 1:
            return False
        if correct_nums == 0 and (correct_count > 0 or wrong_pos_count > 0):
            return False
            
        # Check number size conditions
        if num_size == "both_small" and any(int(guess[i]) >= int(candidate[i]) for i in range(2)):
            return False
        if num_size == "one_small" and sum(1 for i in range(2) if int(guess[i]) < int(candidate[i])) != 1:
            return False
        if num_size == "one_large" and sum(1 for i in range(2) if int(guess[i]) > int(candidate[i])) != 1:
            return False
            
        # Check letters
        if correct_letters == 1 and letter_pos == [2] and guess[2] != candidate[2]:
            return False
        if correct_letters == 1 and letter_pos == "wrong" and guess[3] != candidate[2] and guess[2] != candidate[3]:
            return False
            
    return True

def find_solution():
    # We know 8 is in correct position (from 85BO)
    # We know J is in position 3 (from multiple conditions)
    for n2 in range(10):
        if n2 == 8:  # numbers can't repeat
            continue
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            if letter == 'J':  # letters can't repeat
                continue
            candidate = f"8{n2}J{letter}"
            if validate_against_all_conditions(candidate):
                print(f"Solution found: {candidate}")

find_solution()