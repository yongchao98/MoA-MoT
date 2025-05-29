def check_all_clues(guess):
    n1, n2, l1, l2 = guess
    
    # Clue 1: 79EQ - both numbers too large, both letters incorrect
    if n1 in [7,9] or n2 in [7,9] or l1 == 'E' or l2 == 'Q':
        return False
        
    # Clue 2: 37TK - both numbers incorrect, letters too late
    if n1 == 3 or n2 == 3 or n1 == 7 or n2 == 7:
        return False
    if l1 in 'TK' or l2 in 'TK':
        return False
        
    # Clue 3: 25EO - both numbers incorrect, both letters incorrect
    if n1 in [2,5] or n2 in [2,5] or l1 in 'EO' or l2 in 'EO':
        return False
        
    # Clue 4: 86LZ - one number correct (8) but wrong position, 6 too large
    if n1 == 6 or n2 == 6:
        return False
    if n1 != 8 and n2 != 8:  # 8 must be one of the numbers
        return False
    if l1 in 'LZ' or l2 in 'LZ':
        return False
        
    # Clue 5: 43WQ - both numbers incorrect, letters too late
    if n1 in [4,3] or n2 in [4,3]:
        return False
    if l1 in 'WQ' or l2 in 'WQ':
        return False
        
    # Clue 6: 31AO - both numbers incorrect, A correct but wrong position
    if n1 in [3,1] or n2 in [3,1]:
        return False
    if 'A' not in [l1,l2]:
        return False
        
    # Clue 7: 89GW - numbers too large, G correct and in position
    if n1 == 9 or n2 == 9:
        return False
    if l1 != 'G':  # G must be first letter
        return False
    if l2 == 'W':
        return False
        
    # Basic rules
    if n1 == n2 or l1 == l2:  # No repeats
        return False
        
    return True

# Try all possible combinations
solutions = []
for n1 in range(0,9):  # 0-8
    for n2 in range(0,9):
        for l1 in 'G':  # We know first letter must be G
            for l2 in 'ABCDEFGHIJKLMNOPQRS':  # Excluding letters we know are wrong
                if check_all_clues([n1, n2, l1, l2]):
                    solutions.append([n1, n2, l1, l2])

print(solutions)