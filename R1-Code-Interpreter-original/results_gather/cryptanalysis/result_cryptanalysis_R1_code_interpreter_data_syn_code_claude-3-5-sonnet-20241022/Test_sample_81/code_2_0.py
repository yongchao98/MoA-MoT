from itertools import product

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    num = int(numbers)
    
    # Basic conditions
    if numbers[0] == numbers[1]: return False  # numbers can't repeat
    if numbers[1] != '6': return False  # second number must be 6
    if letters[0] != 'X': return False  # first letter must be X
    
    # Number conditions
    if num <= 2: return False  # too small (from guess 5)
    if num >= 97: return False  # too large (from guess 1)
    
    # From guess 40GU: both numbers incorrect
    if numbers[0] == '4': return False
    
    # From guess 54EW: both numbers incorrect
    if numbers[0] == '5': return False
    
    # From guess 70GF: both numbers incorrect
    if numbers[0] == '7': return False
    
    # From guess 02SV and 02XS: numbers too small
    if num < 30: return False
    
    # From guess 96CL: one number correct (6) and one too large
    if num > 60: return False
    
    # Letter conditions
    wrong_letters = set('YSGFPCLEWNR')  # R is too late (from guess 11)
    if letters[1] in wrong_letters: return False
    
    # From various guesses: letters must be between G and Q
    if not ('G' <= letters[1] <= 'Q'): return False
    
    # From guess 70GF and others: some letters are too early in alphabet
    if letters[1] < 'H': return False
    
    return True

# Generate all possible combinations
valid_combinations = []
for n1 in '3':  # Based on all conditions, first number must be 3
    n2 = '6'    # We know second number is 6
    for l2 in 'HIJKLMOPQ':  # Possible second letters
        guess = [n1, n2, 'X', l2]
        if check_conditions(''.join(guess)):
            valid_combinations.append(guess)

print(valid_combinations)