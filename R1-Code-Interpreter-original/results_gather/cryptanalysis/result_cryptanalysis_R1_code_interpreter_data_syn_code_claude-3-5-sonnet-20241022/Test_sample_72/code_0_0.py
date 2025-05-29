def check_conditions():
    # All possible digits and letters
    digits = list(range(10))
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    # Store potential solutions
    solutions = []
    
    # Try all possible combinations
    for d1 in digits:
        for d2 in digits:
            if d1 == d2:  # No repeating numbers
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:  # No repeating letters
                        continue
                    
                    # Current guess
                    guess = [str(d1), str(d2), l1, l2]
                    
                    # Check all conditions
                    if (
                        # Condition 1: 25WK
                        (d1 != 2 and d2 != 5) and 
                        ((l1 == 'W' or l2 == 'K') != (l1 == 'K' or l2 == 'W')) and
                        
                        # Condition 2: 98DR - numbers too large
                        d1 < 9 and d2 < 8 and
                        
                        # Condition 5: 78LI - one number correct in position
                        ((d1 == 7 and d2 != 8) or (d1 != 7 and d2 == 8)) and
                        
                        # Condition 6: 05KA - one number correct but wrong position
                        ((d1 == 5 and d2 != 0) or (d1 == 0 and d2 != 5) or 
                         (d2 == 5 and d1 != 0) or (d2 == 0 and d1 != 5)) and
                        
                        # Condition 7: 29WH
                        (d1 != 2 and d2 != 9) and
                        ((l1 == 'W' or l2 == 'H') != (l1 == 'H' or l2 == 'W')) and
                        
                        # Condition 9: 72AY - one number correct in position
                        ((d1 == 7 and d2 != 2) or (d1 != 7 and d2 == 2)) and
                        
                        # Condition 10: 45VW - both letters correct in position
                        l1 == 'V' and l2 == 'W'
                    ):
                        solutions.append([str(d1), str(d2), l1, l2])
    
    return solutions

solutions = check_conditions()
print(solutions)