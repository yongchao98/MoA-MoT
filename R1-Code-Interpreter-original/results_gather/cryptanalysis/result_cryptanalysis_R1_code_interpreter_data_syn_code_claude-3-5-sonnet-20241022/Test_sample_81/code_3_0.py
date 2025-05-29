def check_guess(candidate, test_guess, feedback):
    """
    Check if candidate matches the feedback from a test guess
    Returns False if candidate contradicts the feedback
    """
    matches = sum(1 for i in range(4) if candidate[i] == test_guess[i])
    
    # Check number feedback
    if feedback.get('numbers') == 'both_incorrect':
        if candidate[0] == test_guess[0] or candidate[1] == test_guess[1]:
            return False
    elif feedback.get('numbers') == 'one_correct_position':
        if matches != 1:
            return False
            
    # Check letter feedback
    if feedback.get('letters') == 'one_correct_position':
        letter_matches = sum(1 for i in [2,3] if candidate[i] == test_guess[i])
        if letter_matches != 1:
            return False
    elif feedback.get('letters') == 'both_incorrect':
        if candidate[2] == test_guess[2] or candidate[3] == test_guess[3]:
            return False
            
    return True

def check_all_conditions(guess):
    # Basic structure checks
    if guess[1] != '6':  # We know second digit must be 6
        return False
    if guess[2] != 'X':  # We know first letter must be X
        return False
    if int(guess[0] + guess[1]) <= 30 or int(guess[0] + guess[1]) >= 60:
        return False
    
    # Check against all known guesses
    guesses_feedback = [
        ('97YS', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('70GF', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('17YC', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('72XP', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('02SV', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('96CL', {'numbers': 'one_correct_position', 'letters': 'both_incorrect'}),
        ('54EW', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('02XS', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('01GF', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('40GU', {'numbers': 'both_incorrect', 'letters': 'both_incorrect'}),
        ('09XR', {'numbers': 'both_incorrect', 'letters': 'one_correct_position'})
    ]
    
    # Check against each guess
    for test_guess, feedback in guesses_feedback:
        if not check_guess(guess, test_guess, feedback):
            return False
    
    # Additional letter constraints
    if guess[3] in 'YSGFPCLEWNR':  # Known incorrect letters
        return False
    if not 'H' <= guess[3] <= 'Q':  # Must be between H and Q
        return False
        
    return True

# Generate and test all possible combinations
valid_combinations = []
for n1 in '3':  # First number must be 3 based on all constraints
    n2 = '6'    # Second number must be 6
    for l2 in 'HIJKMOQ':  # Possible second letters
        guess = [n1, n2, 'X', l2]
        if check_all_conditions(''.join(guess)):
            valid_combinations.append(guess)

print(valid_combinations)