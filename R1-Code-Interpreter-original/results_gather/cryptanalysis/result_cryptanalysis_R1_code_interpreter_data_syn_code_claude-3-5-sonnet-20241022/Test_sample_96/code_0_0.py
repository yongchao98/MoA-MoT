def check_conditions():
    # All possible digits and letters
    digits = list(range(10))
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    # Store valid combinations
    valid_combinations = []
    
    # Based on condition 7: One number is 3 and it's in correct position (first position)
    # Based on conditions 8,9,10,11,12,13,14: Second number must be less than 4
    
    for d1 in [3]:  # First digit must be 3
        for d2 in range(4):  # Second digit must be less than 4
            if d1 == d2:
                continue
            
            # Based on condition 8: I or J must be in the password but in wrong position
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    
                    combination = [d1, d2, l1, l2]
                    valid = True
                    
                    # Check specific conditions
                    
                    # Condition 7
                    if not (d1 == 3):
                        valid = False
                        
                    # Condition 8
                    if not ((l1 in ['I', 'J'] or l2 in ['I', 'J'])):
                        valid = False
                        
                    # Condition 12
                    if not (d1 == 6 or d2 == 6) and (l1 in ['A', 'B', 'C'] or l2 in ['A', 'B', 'C']):
                        valid = False
                        
                    # Condition 13
                    if not (d1 == 9 or d2 == 9):
                        valid = False
                    
                    # Letters must not be late in alphabet (from condition 1)
                    if ord(l1) > ord('P') or ord(l2) > ord('P'):
                        valid = False
                    
                    if valid:
                        valid_combinations.append([d1, d2, l1, l2])
    
    return valid_combinations

results = check_conditions()
print("Possible combinations:")
for combo in results:
    print(combo)