def is_valid_combination(numbers, letters, combination):
    # Check each guess against the feedback
    guesses = [
        ("38TA", (1, 1, 1, 0)),  # (correct number wrong position, incorrect number too large, correct letter correct position, incorrect letter too late)
        ("10GY", (0, 0, 0, 0)),  # (both numbers too small, both letters too late)
        ("19QH", (0, 0, 0, 0)),  # (both numbers incorrect, both letters too late)
        ("12EF", (0, 0, 1, 0)),  # (both numbers too small, one letter correct wrong position, one letter incorrect)
        ("65TL", (1, 1, 0, 0))   # (correct number wrong position, incorrect number too large, both letters too late)
    ]

    for guess, feedback in guesses:
        guess_numbers = [int(guess[0]), int(guess[1])]
        guess_letters = [guess[2], guess[3]]

        # Check numbers
        correct_number_wrong_position = sum(
            (n in guess_numbers and n != combination[i]) for i, n in enumerate(numbers)
        )
        incorrect_number_too_large = sum(
            (n > max(guess_numbers)) for n in numbers
        )

        # Check letters
        correct_letter_correct_position = sum(
            (l == guess_letters[i]) for i, l in enumerate(letters)
        )
        incorrect_letter_too_late = sum(
            (l > max(guess_letters)) for l in letters
        )

        # Validate feedback
        if (correct_number_wrong_position != feedback[0] or
            incorrect_number_too_large != feedback[1] or
            correct_letter_correct_position != feedback[2] or
            incorrect_letter_too_late != feedback[3]):
            return False

    return True

def find_password():
    possible_numbers = [3, 6, 7, 8]  # Based on feedback, these are the only possible numbers
    possible_letters = ['A', 'B', 'C', 'D', 'E']  # Based on feedback, these are the only possible letters

    for n1 in possible_numbers:
        for n2 in possible_numbers:
            if n1 == n2:
                continue
            for l1 in possible_letters:
                for l2 in possible_letters:
                    if l1 == l2:
                        continue
                    combination = [n1, n2, l1, l2]
                    if is_valid_combination([n1, n2], [l1, l2], combination):
                        return combination

password = find_password()
print(f"<<< {password} >>>")