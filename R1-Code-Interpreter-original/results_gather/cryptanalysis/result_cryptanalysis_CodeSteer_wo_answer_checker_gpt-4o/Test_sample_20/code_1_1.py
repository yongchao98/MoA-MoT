def is_valid_combination(numbers, letters, conditions):
    for guess, feedback in conditions:
        num_correct = 0
        num_incorrect_large = 0
        num_incorrect_small = 0
        letter_correct = 0
        letter_incorrect_late = 0

        # Check numbers
        for i in range(2):
            if numbers[i] == guess[i]:
                num_correct += 1
            elif numbers[i] > guess[i]:
                num_incorrect_large += 1
            else:
                num_incorrect_small += 1

        # Check letters
        for i in range(2):
            if letters[i] == guess[i + 2]:
                letter_correct += 1
            elif letters[i] > guess[i + 2]:
                letter_incorrect_late += 1

        # Validate feedback
        if feedback != (num_correct, num_incorrect_large, num_incorrect_small, letter_correct, letter_incorrect_late):
            return False

    return True

def find_password():
    conditions = [
        ((3, 5, 'P', 'Q'), (1, 1, 0, 0, 2)),
        ((2, 4, 'X', 'E'), (1, 0, 1, 1, 1)),
        ((9, 6, 'S', 'Q'), (0, 2, 0, 0, 2)),
        ((2, 5, 'X', 'C'), (0, 0, 0, 1, 1))
    ]

    # Possible numbers and letters
    possible_numbers = range(10)
    possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

    # Backtracking to find the correct combination
    for n1 in possible_numbers:
        for n2 in possible_numbers:
            if n1 == n2:
                continue
            for l1 in possible_letters:
                for l2 in possible_letters:
                    if l1 == l2:
                        continue
                    if is_valid_combination((n1, n2), (l1, l2), conditions):
                        return [n1, n2, l1, l2]

# Find and print the password
password = find_password()
print(f"<<< {password} >>>")