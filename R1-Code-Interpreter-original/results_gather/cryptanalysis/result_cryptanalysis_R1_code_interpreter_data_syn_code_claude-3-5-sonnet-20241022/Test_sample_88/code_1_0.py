def check_detailed_conditions(guess):
    # Convert guess to string format
    n1, n2, l1, l2 = guess
    
    # Clue 1: 34DU - both numbers and letters wrong
    if n1 in '34' or n2 in '34' or l1 in 'DU' or l2 in 'DU':
        return False
    
    # Clue 2: 24IT - one letter correct but wrong position
    # Since we know V is correct in last position, T must be wrong
    if l1 == 'T' or l2 == 'T' or n1 == '2' or n2 == '2':
        return False
    
    # Clue 5: 07BS - one number correct but wrong position
    # We know 0 is correct in first position, so 7 must be wrong
    if n2 == '7':
        return False
    
    # Clue 6: 42MA - both numbers wrong, letters too early
    if n1 == '4' or n2 == '4' or l1 in 'MA':
        return False
    
    # Clue 8,11: 90,90 - 0 correct in position, 9 too large
    if n1 != '0' or n2 == '9':
        return False
    
    # Clue 12: 10 - 0 correct in position, 1 incorrect
    if n2 == '1':
        return False
    
    # Clue 16: 20KV - 0 correct, V correct, K wrong
    if l2 != 'V' or l1 == 'K':
        return False
    
    # From multiple clues about letters being too early
    if l1 < 'U':
        return False
    
    # From clues about numbers being too large
    if int(n2) >= 8:
        return False
    
    # Clue 7: 38QT - one letter correct but wrong position
    # Since V is correct in last position, Q must be wrong
    if l1 == 'Q':
        return False
    
    # Additional number constraints from various clues
    if n2 in '3489':
        return False
    
    return True

valid_guesses = []
for n2 in '56':  # Remaining possible digits after all constraints
    for l1 in 'WXYZ':  # Remaining possible letters after all constraints
        guess = ['0', n2, l1, 'V']
        if check_detailed_conditions(guess):
            valid_guesses.append(guess)

print(valid_guesses)