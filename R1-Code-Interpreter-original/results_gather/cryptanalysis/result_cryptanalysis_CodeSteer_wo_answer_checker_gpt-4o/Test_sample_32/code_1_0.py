def is_valid_combination(combination):
    # Define the feedback conditions
    conditions = [
        ("19HE", (0, 0, 1, 0)),  # both numbers incorrect; one letter correct and in position; one letter incorrect and too late
        ("75IS", (0, 0, 0, 0)),  # both numbers incorrect; both letters incorrect and too late
        ("58JW", (0, 1, 0, 0)),  # one number correct but wrong position; one number incorrect; both letters incorrect and too late
        ("64XZ", (0, 0, 0, 0)),  # both numbers incorrect; both letters incorrect and too late
        ("07HM", (0, 1, 0, 0)),  # one number correct but wrong position; one number incorrect; both letters incorrect and too late
        ("35WS", (0, 0, 0, 0)),  # both numbers incorrect; both letters incorrect and too late
        ("47HL", (0, 0, 0, 0)),  # both numbers incorrect; both letters incorrect and too late
        ("23QE", (0, 0, 1, 0)),  # both numbers incorrect; one letter correct and in position; one letter incorrect and too late
        ("48AJ", (0, 1, 1, 0)),  # one number correct but wrong position; one number incorrect; one letter correct and in position; one letter incorrect and too late
    ]

    for guess, feedback in conditions:
        num_correct_pos = sum(combination[i] == guess[i] for i in range(2))
        num_correct_wrong_pos = sum(combination[i] in guess[:2] and combination[i] != guess[i] for i in range(2))
        letter_correct_pos = sum(combination[i+2] == guess[i+2] for i in range(2))
        letter_correct_wrong_pos = sum(combination[i+2] in guess[2:] and combination[i+2] != guess[i+2] for i in range(2))

        # Debug statement to show the comparison of feedback
        print(f"Checking combination {combination} against guess {guess}:")
        print(f"Expected feedback: {feedback}")
        print(f"Calculated feedback: {(num_correct_pos, num_correct_wrong_pos, letter_correct_pos, letter_correct_wrong_pos)}")

        if (num_correct_pos, num_correct_wrong_pos, letter_correct_pos, letter_correct_wrong_pos) != feedback:
            return False
    return True

def find_password():
    numbers = ['0', '8']
    letters = ['Q', 'A']
    for num1 in numbers:
        for num2 in numbers:
            if num1 == num2:
                continue
            for letter1 in letters:
                for letter2 in letters:
                    if letter1 == letter2:
                        continue
                    combination = [num1, num2, letter1, letter2]
                    print(f"Trying combination: {combination}")  # Debug statement
                    if is_valid_combination(combination):
                        print(f"Valid combination found: {combination}")  # Debug statement
                        return combination

password = find_password()
print(f"<<< {password} >>>")