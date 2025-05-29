def check_number_conditions(n1, n2):
    # Numbers can't repeat
    if n1 == n2:
        return False
        
    # Second number must be 6 (from guess 96CL where one number is correct in position)
    if n2 != '6':
        return False
        
    num = int(n1 + n2)
    # Must be > 02 (too small) but < 97 (too large)
    if num <= 30 or num >= 60:  # Refined range based on all guesses
        return False
        
    # Can't be any of these first digits (from guesses where both numbers were wrong)
    if n1 in ['0', '1', '2', '4', '5', '7', '9']:
        return False
        
    return True

def check_letter_conditions(l1, l2):
    # First letter must be X (from guess 09XR where X is correct and R is too late)
    if l1 != 'X':
        return False
        
    # Letters that are known to be incorrect from various guesses
    wrong_letters = set('YSGFPCLEWNR')
    if l2 in wrong_letters:
        return False
        
    # Must be between H and Q (from various feedbacks about too early/late)
    if not ('H' <= l2 <= 'Q'):
        return False
        
    # Additional constraints from guesses where both letters were wrong
    if l2 in ['H', 'I', 'J', 'K', 'O', 'Q']:  # These would make some guesses partially correct
        return False
        
    return True

# Generate all possible combinations
valid_combinations = []
for n1 in '38':  # Only possible first digits based on all constraints
    for n2 in '6':  # We know second digit is 6
        if check_number_conditions(n1, n2):
            for l2 in 'M':  # After all constraints, only M remains possible
                if check_letter_conditions('X', l2):
                    valid_combinations.append([n1, n2, 'X', l2])

print(valid_combinations)