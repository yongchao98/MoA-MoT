def is_valid_combination(combination):
    # Debugging: Print the current combination being checked
    print(f"Checking combination: {combination}")

    # Condition 1: 98SO
    # Both numbers are incorrect and too large; one letter is correct and in the correct position; one letter is incorrect.
    if not (combination[0] < 9 and combination[1] < 8 and
            ((combination[2] == 'S' and combination[3] != 'O') or
             (combination[2] != 'S' and combination[3] == 'O'))):
        return False

    # Condition 2: 34ZW
    # One number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect and too late in the alphabet.
    if not ((combination[0] == 3 or combination[1] == 3) and
            (combination[0] != 4 and combination[1] != 4) and
            combination[2] < 'Z' and combination[3] < 'W'):
        return False

    # Condition 3: 91HL
    # Both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet.
    if not ((combination[2] == 'H' and combination[3] != 'L') or
            (combination[2] != 'H' and combination[3] == 'L')):
        return False

    # Condition 4: 47OB
    # One number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect.
    if not ((combination[0] == 4 and combination[1] != 7) or
            (combination[0] != 4 and combination[1] == 7)):
        return False

    # Condition 5: 65NA
    # One number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect.
    if not ((combination[0] == 6 and combination[1] != 5) or
            (combination[0] != 6 and combination[1] == 5)):
        return False

    return True

def find_password():
    numbers = [2, 3, 4, 5, 6, 7]
    letters = ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'S']

    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    combination = [n1, n2, l1, l2]
                    if is_valid_combination(combination):
                        print(f"Valid combination found: {combination}")
                        return combination

password = find_password()
print(f"<<< {password} >>>")