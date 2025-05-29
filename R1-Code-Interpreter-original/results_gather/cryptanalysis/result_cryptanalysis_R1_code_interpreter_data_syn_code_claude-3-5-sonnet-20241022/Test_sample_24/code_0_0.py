from itertools import permutations

def check_guess(guess, actual):
    # Returns feedback matching the format in the clues
    correct_pos = 0
    correct_wrong_pos = 0
    for i in range(4):
        if guess[i] == actual[i]:
            correct_pos += 1
        elif guess[i] in actual:
            correct_wrong_pos += 1
    return correct_pos, correct_wrong_pos

def is_letter_too_early(guess_letter, actual_letter):
    return guess_letter < actual_letter

def validate_password(password):
    # Check against all three clues
    
    # Clue 1: 95WC
    # - one number correct wrong pos
    # - one number incorrect
    # - one letter correct pos
    # - one letter incorrect and too early
    clue1 = "95WC"
    correct_pos1, correct_wrong_pos1 = check_guess(clue1, password)
    if not (correct_pos1 == 1 and correct_wrong_pos1 == 1):
        return False
    
    # Clue 2: 12XN
    # - one number correct pos
    # - one number incorrect
    # - both letters incorrect
    clue2 = "12XN"
    correct_pos2, correct_wrong_pos2 = check_guess(clue2, password)
    if not (correct_pos2 == 1 and correct_wrong_pos2 == 0):
        return False
    
    # Clue 3: 92LS
    # - one number correct wrong pos
    # - one number incorrect
    # - one letter correct wrong pos
    # - one letter incorrect
    clue3 = "92LS"
    correct_pos3, correct_wrong_pos3 = check_guess(clue3, password)
    if not (correct_pos3 == 0 and correct_wrong_pos3 == 2):
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