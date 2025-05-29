def verify_combination(combo):
    # Convert combo to list of strings
    combo = [str(x) for x in combo]
    
    # Check against guess 1: 64RY - all incorrect
    g1 = ['6', '4', 'R', 'Y']
    if any(x in combo for x in g1):
        return False
    
    # Check against guess 2: 51MB - all incorrect
    g2 = ['5', '1', 'M', 'B']
    if any(x in combo for x in g2):
        return False
    
    # Check against guess 3: 29PV - all incorrect
    g3 = ['2', '9', 'P', 'V']
    if any(x in combo for x in g3):
        return False
    
    # Check against guess 4: 71EW
    g4 = ['7', '1', 'E', 'W']
    
    # One number correct but wrong position (7)
    # Both letters correct but wrong positions
    correct_num = sum(1 for i in range(2) if combo[i] == '7')
    if correct_num != 1:
        return False
    if combo[0] == '7':  # 7 must be in different position than guess 4
        return False
        
    # Check letters (E,W must be present but in different positions than guess 4)
    if not ('E' in combo and 'W' in combo):
        return False
    if combo[2] == 'E' and combo[3] == 'W':  # Can't be in same positions as guess 4
        return False
        
    return True

# Test all possible combinations
candidates = [['0', '7', 'W', 'E'], ['3', '7', 'W', 'E'], ['8', '7', 'W', 'E']]
for combo in candidates:
    if verify_combination(combo):
        print(f"Valid combination found: {combo}")