def check_conditions(guess):
    n1, n2, l1, l2 = guess
    
    # Helper function to convert string numbers to integers
    n1, n2 = int(n1), int(n2)
    
    # Condition 1: 93YW
    # Both numbers correct and in position, letters too late
    if not (n1 == 9 and n2 == 3):
        return False
    if not (l1 < 'Y' and l2 < 'Y'):
        return False
    
    # Condition 5: 46GK
    # Numbers incorrect, one letter wrong pos, one incorrect
    if 'G' not in [l1, l2] or 'K' not in [l1, l2]:
        return False
    
    # Condition 10: 43KU
    # One number correct pos (3), one letter wrong pos
    if 'K' in [l1, l2]:
        if l2 == 'K':  # K must be in wrong position
            return False
    
    # Additional checks from other conditions
    # Letters can't be too early (from condition 9: 61BE)
    if l1 <= 'E' or l2 <= 'E':
        return False
    
    # Letters can't be certain values from various conditions
    forbidden_letters = {'Y', 'W', 'H', 'F', 'Z', 'R', 'P', 'V', 'N', 'C', 'B', 'E', 'U'}
    if l1 in forbidden_letters or l2 in forbidden_letters:
        return False
    
    return True

# Generate and test all possible combinations
numbers = ['9', '3']  # We know these from condition 1
letters = [chr(i) for i in range(ord('A'), ord('Z')+1)]
valid_combinations = []

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # Letters can't repeat
            guess = [numbers[0], numbers[1], l1, l2]
            if check_conditions(guess):
                valid_combinations.append(guess)

print(valid_combinations)