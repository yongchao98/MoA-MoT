from itertools import product
import string

def is_valid_guess(guess, numbers, letters):
    # Check if numbers and letters match all conditions from the clues
    
    # From clue 1 (34DU): both numbers wrong, both letters wrong
    if guess[0] in ['3', '4'] or guess[1] in ['3', '4'] or guess[2] in 'DU' or guess[3] in 'DU':
        return False
    
    # From clue 8,11,12,14,16: first digit is 0
    if guess[0] != '0':
        return False
    
    # From clue 16: V is correct in position 4
    if guess[3] != 'V':
        return False
    
    # From clues about numbers being too large
    if int(guess[1]) >= 8:  # Based on clues about numbers being too large
        return False
        
    # From clue 16: one number correct and in position, other incorrect
    if guess[0] != '0':
        return False
        
    # From various clues about letters being too early
    early_letters = set('ABCDEFGHIJKLMNOPQRST')
    if guess[2] in early_letters:
        return False
    
    return True

# Generate all possible combinations
digits = '0123456789'
letters = string.ascii_uppercase

valid_guesses = []
for n1, n2 in product(digits, digits):
    if n1 == n2:  # Skip if numbers repeat
        continue
    for l1, l2 in product(letters, letters):
        if l1 == l2:  # Skip if letters repeat
            continue
        guess = [n1, n2, l1, l2]
        if is_valid_guess(guess, [n1, n2], [l1, l2]):
            valid_guesses.append(guess)

print(valid_guesses)