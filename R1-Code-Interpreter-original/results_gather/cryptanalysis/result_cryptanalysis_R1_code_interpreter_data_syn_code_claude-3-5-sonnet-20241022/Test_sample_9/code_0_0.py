def find_password():
    # All possible numbers (excluding known incorrect ones)
    possible_numbers = [0,3,7,8]  # Excluding 1,2,4,5,6,9
    
    # We know E,W are the letters, just in wrong positions
    letters = ['E', 'W']
    
    # From guess 4 (71EW), we know 7 or 1 is correct but in wrong position
    # Since 1 was proven wrong in guess 2, 7 must be the correct number
    
    valid_combinations = []
    
    # Try all possible combinations
    for n1 in possible_numbers:
        for n2 in possible_numbers:
            if n1 != n2:  # Numbers can't repeat
                for l1 in letters:
                    for l2 in letters:
                        if l1 != l2:  # Letters can't repeat
                            candidate = [str(n1), str(n2), l1, l2]
                            
                            # Check against all conditions
                            # Condition from guess 4: 7 must be in different position than in 71EW
                            if '7' in [str(n1), str(n2)]:
                                if str(n1) == '7':  # 7 was first in guess 4, so it must be second here
                                    continue
                                
                                # Check if E,W are in different positions than in guess 4
                                if l1 == 'E' and l2 == 'W':  # Same positions as guess 4
                                    continue
                                    
                                valid_combinations.append(candidate)
    
    # Print all valid combinations
    for combo in valid_combinations:
        print(combo)

find_password()