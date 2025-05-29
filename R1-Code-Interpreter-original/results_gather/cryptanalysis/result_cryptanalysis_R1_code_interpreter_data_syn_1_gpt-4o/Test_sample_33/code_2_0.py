# Possible numbers and letters based on deductions
numbers = [8, 9]
letters = ['Y', 'Q']

# Check each condition
def check_conditions():
    # Condition 1: 35WJ
    if not (numbers[0] == 3 or numbers[1] == 5):
        return False
    if not (letters[0] != 'W' and letters[1] != 'J'):
        return False

    # Condition 2: 67VX
    if 6 in numbers or 7 in numbers:
        return False
    if 'V' in letters or 'X' in letters:
        return False

    # Condition 3: 37QK
    if 3 in numbers or 7 in numbers:
        return False
    if not (letters[0] == 'Q' and letters[1] != 'K'):
        return False

    # Condition 4: 75GT
    if not (numbers[0] == 7 or numbers[1] == 5):
        return False
    if 'G' in letters or 'T' in letters:
        return False

    # Condition 5: 34QE
    if 3 in numbers or 4 in numbers:
        return False
    if 'Q' in letters or 'E' in letters:
        return False

    # Condition 6: 27CU
    if 2 in numbers or 7 in numbers:
        return False
    if 'C' in letters or 'U' in letters:
        return False

    # Condition 7: 56FA
    if not (numbers[0] == 5 or numbers[1] == 6):
        return False
    if 'F' in letters or 'A' in letters:
        return False

    # Condition 8: 84YT
    if not (numbers[0] == 8 and numbers[1] != 4):
        return False
    if not (letters[0] == 'Y' and letters[1] != 'T'):
        return False

    return True

# Check if the conditions are satisfied
if check_conditions():
    print([numbers[0], numbers[1], letters[0], letters[1]])
else:
    print("No valid combination found.")