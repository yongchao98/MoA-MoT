def verify_guess(guess, numbers, letters):
    # Returns True if the guess matches all conditions
    d1, d2, l1, l2 = guess
    
    # Condition 1: 16QS - both numbers incorrect, letters too late
    if d1 in [1, 6] or d2 in [1, 6]:
        return False
    
    # Condition 2: 23YQ - one number correct wrong position
    if d1 in [2, 3] and d2 in [2, 3]:
        return False
        
    # Condition 7: 37HS - first number correct position, second too large
    if d1 != 3 or d2 >= 7:
        return False
        
    # Condition 8: 58IJ - numbers too large, one letter correct wrong position
    if d1 >= 5 or d2 >= 5:
        return False
    if not ('I' in letters or 'J' in letters):
        return False
        
    # Condition 11: 43PR - one number correct wrong position
    if d1 == 4 or d2 == 4:
        return False
        
    # Condition 12: 60CA - one number correct position
    if l1 in ['A', 'B', 'C'] or l2 in ['A', 'B', 'C']:
        return False
        
    # Letters can't be too late in alphabet
    if ord(l1) > ord('P') or ord(l2) > ord('P'):
        return False
        
    return True

# Generate all possible combinations
valid_combinations = []
for d2 in range(4):  # Second digit must be less than 4
    if d2 == 3:  # Can't be 3 as it's already first digit
        continue
    for l1 in [chr(i) for i in range(ord('D'), ord('P')+1)]:  # D-P
        for l2 in [chr(i) for i in range(ord('D'), ord('P')+1)]:
            if l1 == l2:
                continue
            if l1 != 'I' and l2 != 'I':  # One letter must be I
                continue
            guess = [3, d2, l1, l2]
            if verify_guess(guess, [3, d2], [l1, l2]):
                valid_combinations.append(guess)

print("Final possible combinations:")
for combo in valid_combinations:
    print(combo)