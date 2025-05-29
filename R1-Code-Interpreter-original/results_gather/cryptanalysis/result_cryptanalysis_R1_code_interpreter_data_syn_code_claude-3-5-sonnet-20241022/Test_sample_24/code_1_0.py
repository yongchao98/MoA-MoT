from itertools import permutations

def check_guess(guess, actual):
    # Returns positions of correct and wrong position matches
    correct_pos = []
    correct_wrong_pos = []
    for i in range(4):
        if guess[i] == actual[i]:
            correct_pos.append(i)
        elif guess[i] in actual:
            correct_wrong_pos.append(i)
    return correct_pos, correct_wrong_pos

def is_letter_too_early(letter1, letter2):
    return ord(letter1) < ord(letter2)

def validate_password(password):
    # Clue 1: 95WC
    clue1 = "95WC"
    correct_pos1, correct_wrong_pos1 = check_guess(clue1, password)
    # Must have exactly one letter correct in position (W or C)
    # One number correct but wrong position (9 or 5)
    # One letter must be too early in alphabet
    letter_correct = False
    number_wrong_pos = False
    too_early_found = False
    
    if len(correct_pos1) != 1 or len(correct_wrong_pos1) != 1:
        return False
    
    # Check if the correct position is a letter
    if correct_pos1[0] >= 2:
        letter_correct = True
    
    # Check if wrong position is a number
    if correct_wrong_pos1[0] < 2:
        number_wrong_pos = True
    
    # For the incorrect letter, check if it's too early
    actual_letters = [password[2], password[3]]
    guess_letters = [clue1[2], clue1[3]]
    for g in guess_letters:
        if g not in password and any(is_letter_too_early(g, a) for a in actual_letters):
            too_early_found = True
    
    if not (letter_correct and number_wrong_pos and too_early_found):
        return False

    # Clue 2: 12XN
    clue2 = "12XN"
    correct_pos2, correct_wrong_pos2 = check_guess(clue2, password)
    # Must have exactly one number correct in position
    # Both letters must be incorrect
    if len(correct_pos2) != 1 or len(correct_wrong_pos2) != 0:
        return False
    if correct_pos2[0] >= 2:  # The correct position must be a number
        return False
    if clue2[2] in password[2:] or clue2[3] in password[2:]:  # Letters must be incorrect
        return False

    # Clue 3: 92LS
    clue3 = "92LS"
    correct_pos3, correct_wrong_pos3 = check_guess(clue3, password)
    # One number correct but wrong position
    # One letter correct but wrong position
    if len(correct_pos3) != 0 or len(correct_wrong_pos3) != 2:
        return False
    
    # Must have one number and one letter in wrong positions
    number_wrong = False
    letter_wrong = False
    for pos in correct_wrong_pos3:
        if pos < 2:
            number_wrong = True
        else:
            letter_wrong = True
    
    if not (number_wrong and letter_wrong):
        return False

    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

valid_passwords = []

# Try all possible combinations
for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        password = n1 + n2 + l1 + l2
        if validate_password(password):
            valid_passwords.append([n1, n2, l1, l2])

print(valid_passwords)