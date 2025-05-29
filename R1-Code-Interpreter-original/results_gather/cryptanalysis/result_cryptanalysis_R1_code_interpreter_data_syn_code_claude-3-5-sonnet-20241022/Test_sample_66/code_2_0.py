import string

def check_guess(numbers, letters):
    # Numbers must be 1,8
    if numbers != [1, 8]:
        return False
    
    # Second letter must be E (from guess 6)
    if letters[1] != 'E':
        return False
    
    # First letter must be N (from guess 8 and position deduction)
    if letters[0] != 'N':
        return False
    
    # Letters must be before 'S' in alphabet
    for letter in letters:
        if letter >= 'S':
            return False
    
    # Verify against guess 8 (87NM)
    # N must be correct but in wrong position when in second position
    if letters == ['M', 'N']:
        return False
    
    # Verify against guess 11 (18DN)
    # Both letters must be incorrect in these positions
    if letters[0] == 'D' or letters[1] == 'N':
        return False
    
    return True

# Generate all possible combinations
valid_numbers = [1, 8]
valid_combinations = []

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