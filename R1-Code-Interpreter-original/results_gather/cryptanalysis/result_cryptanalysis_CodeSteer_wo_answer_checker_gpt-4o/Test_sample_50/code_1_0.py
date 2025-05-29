def is_valid_combination(combination):
    feedback_rules = [
        ("15YA", (1, 1, 1, 0)),  # one number correct but wrong position, one number too large, one letter correct and correct position, one letter incorrect and too early
        ("07JZ", (0, 0, 0, 0)),  # both numbers and letters incorrect
        ("70RC", (0, 0, 0, 0)),  # both numbers incorrect, both letters too early
        ("70WL", (0, 0, 0, 0)),  # both numbers and letters incorrect
        ("54XP", (0, 0, 0, 0)),  # both numbers too large, both letters incorrect
        ("59IT", (0, 0, 1, 1)),  # both numbers too large, one letter correct and correct position, one letter too early
        ("41WR", (1, 0, 0, 0)),  # one number correct and correct position, one number too large, both letters incorrect
        ("60RF", (0, 0, 0, 0)),  # both numbers incorrect, both letters too early
        ("63OA", (0, 0, 0, 0)),  # both numbers too large, both letters too early
    ]

    for guess, feedback in feedback_rules:
        guess_numbers = guess[:2]
        guess_letters = guess[2:]

        # Check numbers
        correct_number_position = sum(1 for i in range(2) if combination[i] == guess_numbers[i])
        correct_number_wrong_position = sum(1 for i in range(2) if combination[i] in guess_numbers and combination[i] != guess_numbers[i])
        
        if correct_number_position != feedback[0] or correct_number_wrong_position != feedback[1]:
            return False

        # Check letters
        correct_letter_position = sum(1 for i in range(2, 4) if combination[i] == guess_letters[i - 2])
        correct_letter_wrong_position = sum(1 for i in range(2, 4) if combination[i] in guess_letters and combination[i] != guess_letters[i - 2])
        
        if correct_letter_position != feedback[2] or correct_letter_wrong_position != feedback[3]:
            return False

    return True

def find_password():
    numbers = "0123456789"
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    combination = n1 + n2 + l1 + l2
                    if is_valid_combination(combination):
                        return [n1, n2, l1, l2]

password = find_password()
print(f"<<< {password} >>>")