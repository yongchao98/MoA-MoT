import string

def check_guess(numbers, letters):
    # Known conditions
    if numbers != [1, 8]:
        return False
    
    # Letters must be before 'S' in alphabet
    for letter in letters:
        if letter > 'S':
            return False
            
    # Cannot be any of the known incorrect letters
    wrong_letters = set('FYKSQBSUYWNMDL')
    if any(l in wrong_letters for l in letters):
        return False
    
    # From guess 6: one of U,E is correct in position
    if letters[1] != 'E':  # We know E must be in second position because U is too late
        return False
        
    # From guess 8: N must be present but not in final position
    if letters[0] != 'N':  # Must be in first position as we know E is in second
        return False
        
    return True

# Generate all possible combinations
valid_numbers = [1, 8]
valid_letters = []

for l1 in string.ascii_uppercase:
    for l2 in string.ascii_uppercase:
        if l1 != l2:  # Letters can't repeat
            if check_guess(valid_numbers, [l1, l2]):
                valid_letters.append([l1, l2])

print(f"Valid combination: {valid_numbers + valid_letters[0]}")