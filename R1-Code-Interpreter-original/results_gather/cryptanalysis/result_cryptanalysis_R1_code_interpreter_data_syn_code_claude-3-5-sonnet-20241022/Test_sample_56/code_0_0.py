import itertools
import string

def check_guess(guess, numbers, letters):
    # Returns True if the guess is consistent with all conditions
    
    # Condition 1: 93YW - both numbers correct and in position, letters too late
    if not (guess[0] == '9' and guess[1] == '3' and 
            all(l1 < l2 for l1, l2 in zip(guess[2:], ['Y', 'W']))):
        return False
    
    # Condition 2: 34PL - one number wrong pos, one incorrect, letters incorrect
    if not ((('3' in guess[:2] and guess.index('3') != 1) or 
             ('4' in guess[:2] and guess.index('4') != 0)) and
            not ('3' in guess[:2] and '4' in guess[:2])):
        return False
    
    # Condition 3: 14YH - all incorrect
    if '1' in guess[:2] or '4' in guess[:2] or 'Y' in guess[2:] or 'H' in guess[2:]:
        return False
    
    # Condition 4: 70FZ - all incorrect
    if '7' in guess[:2] or '0' in guess[:2] or 'F' in guess[2:] or 'Z' in guess[2:]:
        return False
    
    # Condition 5: 46GK - numbers incorrect, one letter wrong pos, one incorrect
    if '4' in guess[:2] or '6' in guess[:2]:
        return False
    if not (('G' in guess[2:] and guess.index('G') != 2) or 
            ('K' in guess[2:] and guess.index('K') != 3)):
        return False
    
    # Condition 6: 32RP - one number wrong pos, one too small, letters incorrect
    if not ('3' in guess[:2] and guess.index('3') != 0):
        return False
    if '2' in guess[:2] or 'R' in guess[2:] or 'P' in guess[2:]:
        return False
    
    # Condition 7: 63VH - one number correct pos, one incorrect
    if not (('6' in guess[:2] and guess.index('6') == guess[:2].index('6')) or 
            ('3' in guess[:2] and guess.index('3') == guess[:2].index('3'))):
        return False
    if 'V' in guess[2:] or 'H' in guess[2:]:
        return False
    
    # Condition 8: 13NC - one number correct pos, one too small
    if not ('3' in guess[:2] and guess.index('3') == guess[:2].index('3')):
        return False
    if '1' in guess[:2] or 'N' in guess[2:] or 'C' in guess[2:]:
        return False
    
    # Condition 9: 61BE - numbers incorrect, letters too early
    if '6' in guess[:2] or '1' in guess[:2]:
        return False
    if 'B' in guess[2:] or 'E' in guess[2:]:
        return False
    
    # Condition 10: 43KU - one number correct pos, one incorrect, one letter wrong pos
    if not ('3' in guess[:2] and guess.index('3') == guess[:2].index('3')):
        return False
    if not (('K' in guess[2:] and guess.index('K') != 2) or 
            ('U' in guess[2:] and guess.index('U') != 3)):
        return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = string.ascii_uppercase

valid_combinations = []
for nums in itertools.permutations(numbers, 2):
    for lets in itertools.permutations(letters, 2):
        guess = list(nums + lets)
        if check_guess(guess, numbers, letters):
            valid_combinations.append(guess)

print(valid_combinations)