def analyze_conditions():
    # Known facts:
    # 1. From condition 8 (45TG): both numbers are 4,5 in correct positions
    # 2. From condition 6 (64LY): L is in correct position
    # 3. Let's analyze each condition systematically
    
    numbers = [4, 5]  # confirmed from condition 8
    
    # Track possible letters for each position
    possible_letters = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    
    # From condition 6 (64LY):
    # - L is correct and in position
    # - Y is too late
    third_position = 'L'  # confirmed
    
    # Remove letters we know are too early or too late
    too_early = set('EFHIBKG')
    too_late = set('PZVWXNY')
    possible_letters -= too_early
    possible_letters -= too_late
    
    # Test remaining combinations
    solutions = []
    for letter in possible_letters:
        if letter != 'L':  # Skip L as it's in third position
            candidate = [4, 5, 'L', letter]
            # Verify against all conditions
            valid = True
            
            # Check against condition 4 (58QD)
            if letter == 'D':
                valid = False
                
            # Check against condition 6 (64LY)
            if letter == 'Y':
                valid = False
                
            # Check against condition 7 (48HI)
            if letter in 'HI':
                valid = False
                
            # Check against condition 8 (45TG)
            if letter in 'TG':
                valid = False
            
            if valid:
                solutions.append(candidate)
    
    print("Possible solutions:")
    for solution in solutions:
        print(solution)

analyze_conditions()