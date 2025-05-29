import string

def check_guess(numbers, letters):
    # Known conditions
    if numbers != [1, 8]:
        return False
    
    # Letters must be before 'S' in alphabet
    for letter in letters:
        if letter >= 'S':
            return False
            
    # Cannot be any of the known incorrect letters
    wrong_letters = set('FYKSQBSUYWNMDL')
    if any(l in wrong_letters for l in letters):
        return False
    
    # Check specific conditions from guesses
    if letters[1] != 'E':  # From guess 6
        return False
        
    if letters[0] == 'N':  # N can't be in first position as it's wrong in guess 8
        return False
        
    return True

# Generate all possible combinations
valid_numbers = [1, 8]
valid_combinations = []

# Try all possible letters before 'S'
for l1 in string.ascii_uppercase:
    if l1 >= 'S':
        continue
    for l2 in string.ascii_uppercase:
        if l2 >= 'S':
            continue
        if l1 != l2:  # Letters can't repeat
            if check_guess(valid_numbers, [l1, l2]):
                valid_combinations.append([1, 8, l1, l2])

print("Valid combinations:")
for combo in valid_combinations:
    print(combo)