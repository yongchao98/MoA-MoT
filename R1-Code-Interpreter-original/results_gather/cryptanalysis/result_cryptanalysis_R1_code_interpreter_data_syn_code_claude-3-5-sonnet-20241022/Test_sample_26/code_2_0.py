def is_valid_combination(n1, n2, l1, l2):
    # Helper function to check letter position feedback
    def check_letter_position(guess_letters, actual_letters, feedback):
        if feedback == "both_wrong":
            return not any(g in actual_letters for g in guess_letters)
        elif feedback == "one_correct_wrong_pos":
            correct_count = sum(1 for g in guess_letters if g in actual_letters)
            different_positions = (guess_letters[0] != actual_letters[0]) and (guess_letters[1] != actual_letters[1])
            return correct_count == 1 and different_positions

    # Helper function to check number position feedback
    def check_number_position(guess_nums, actual_nums, feedback):
        if feedback == "both_wrong":
            return not any(str(g) == str(a) for g, a in zip(guess_nums, actual_nums))
        elif feedback == "one_correct_pos":
            return sum(1 for g, a in zip(guess_nums, actual_nums) if str(g) == str(a)) == 1
        elif feedback == "one_correct_wrong_pos":
            return sum(1 for g in guess_nums if str(g) in [str(n1), str(n2)]) == 1 and \
                   not any(str(g) == str(a) for g, a in zip(guess_nums, actual_nums))

    # Convert numbers to strings for easier comparison
    nums = [str(n1), str(n2)]
    letters = [l1, l2]

    # Check Condition 1: 83VM
    if not (check_number_position(['8', '3'], nums, "both_wrong") and
            check_letter_position(['V', 'M'], letters, "one_correct_wrong_pos")):
        return False

    # Check Condition 2: 40WL
    if not (check_number_position(['4', '0'], nums, "both_wrong") and
            check_letter_position(['W', 'L'], letters, "both_wrong")):
        return False

    # Check Condition 3: 50VD
    if not (check_number_position(['5', '0'], nums, "both_wrong") and
            check_letter_position(['V', 'D'], letters, "one_correct_wrong_pos")):
        return False

    # Check Condition 4: 40MC
    if not (check_number_position(['4', '0'], nums, "both_wrong") and
            check_letter_position(['M', 'C'], letters, "both_wrong")):
        return False

    # Check Condition 5: 18JC
    if not (check_number_position(['1', '8'], nums, "one_correct_pos") and
            check_letter_position(['J', 'C'], letters, "both_wrong")):
        return False
    if n1 != 1:  # First number must be 1
        return False

    # Check Condition 6: 41TW
    if not (check_number_position(['4', '1'], nums, "one_correct_wrong_pos") and
            check_letter_position(['T', 'W'], letters, "both_wrong")):
        return False

    # Check Condition 7: 07BF
    if not (check_number_position(['0', '7'], nums, "both_wrong") and
            check_letter_position(['B', 'F'], letters, "one_correct_wrong_pos")):
        return False

    return True

# Test all possible combinations
solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2:
                        if is_valid_combination(n1, n2, l1, l2):
                            solutions.append([n1, n2, l1, l2])

print(solutions)