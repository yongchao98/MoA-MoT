def is_valid_combination(numbers, letters, guess, feedback):
    num_correct_pos = 0
    num_correct_wrong_pos = 0
    num_incorrect_too_small = 0
    num_incorrect_too_large = 0
    letter_correct_pos = 0
    letter_incorrect_too_early = 0
    letter_incorrect_too_late = 0

    # Check numbers
    for i in range(2):
        if numbers[i] == guess[i]:
            num_correct_pos += 1
        elif numbers[i] in guess[:2]:
            num_correct_wrong_pos += 1
        elif numbers[i] < guess[i]:
            num_incorrect_too_small += 1
        else:
            num_incorrect_too_large += 1

    # Check letters
    for i in range(2):
        if letters[i] == guess[i + 2]:
            letter_correct_pos += 1
        elif letters[i] < guess[i + 2]:
            letter_incorrect_too_early += 1
        else:
            letter_incorrect_too_late += 1

    # Validate against feedback
    expected_feedback = (num_correct_pos, num_correct_wrong_pos, num_incorrect_too_small, num_incorrect_too_large, letter_correct_pos, letter_incorrect_too_early, letter_incorrect_too_late)
    if feedback == expected_feedback:
        return True
    return False

def find_password():
    feedbacks = [
        ((2, 5, 'Y', 'N'), (0, 1, 1, 0, 0, 0, 2)),
        ((4, 2, 'C', 'Q'), (0, 0, 2, 0, 0, 0, 2)),
        ((7, 6, 'K', 'C'), (1, 0, 0, 1, 1, 1, 0)),
        ((2, 6, 'S', 'G'), (1, 0, 1, 0, 1, 0, 1))
    ]

    for num1 in range(10):
        for num2 in range(10):
            if num1 == num2:
                continue
            for letter1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for letter2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if letter1 == letter2:
                        continue
                    valid = True
                    for guess, feedback in feedbacks:
                        if not is_valid_combination([num1, num2], [letter1, letter2], guess, feedback):
                            valid = False
                            break
                    if valid:
                        return [str(num1), str(num2), letter1, letter2]

password = find_password()
print(f"<<< {password} >>>")