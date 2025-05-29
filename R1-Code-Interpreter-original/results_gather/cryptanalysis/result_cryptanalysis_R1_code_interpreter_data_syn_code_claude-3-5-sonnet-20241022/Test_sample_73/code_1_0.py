def verify_all_conditions(pwd):
    n1, n2, l1, l2 = pwd
    
    conditions = [
        # 1: 74SE - both numbers incorrect, one letter correct position, one too early
        lambda: n1 != '7' and n2 != '4' and ((l1 == 'S') != (l2 == 'E')),
        
        # 2: 76OX - both numbers and letters incorrect
        lambda: n1 != '7' and n2 != '6' and l1 != 'O' and l2 != 'X',
        
        # 3: 94LJ - both numbers incorrect, letters too early
        lambda: n1 != '9' and n2 != '4' and l1 != 'L' and l2 != 'J',
        
        # 4: 28HI - one number correct position (8), one too small, letters too early
        lambda: n2 == '8' and int(n1) > 2 and l1 != 'H' and l2 != 'I',
        
        # 5: 27LJ - both numbers incorrect, letters too early
        lambda: n1 != '2' and n2 != '7' and l1 != 'L' and l2 != 'J',
        
        # 6: 36WF - both numbers and letters incorrect
        lambda: n1 != '3' and n2 != '6' and l1 != 'W' and l2 != 'F',
        
        # 7: 71IS - one letter correct wrong position
        lambda: n1 != '7' and n2 != '1' and ((l2 == 'S' and l1 != 'I') or (l1 != 'S' and l2 != 'I')),
        
        # 8: 12YB - both numbers too small
        lambda: int(n1) > 1 and int(n2) > 2,
        
        # 9: 87WE - one number wrong position
        lambda: ((n1 == '7' or n2 == '8') and not (n1 == '8' and n2 == '7')) and l1 != 'W' and l2 != 'E',
        
        # 10: 25YP - one number wrong position, one too small
        lambda: int(n1) > 2 and n2 != '5' and l1 != 'Y' and l2 != 'P',
        
        # 11: 84UL - one number wrong position, one too small
        lambda: int(n1) > 4 and n2 != '4' and l1 != 'U' and l2 != 'L',
        
        # 12: 30VF - both numbers too small
        lambda: int(n1) > 3 and int(n2) > 0,
        
        # 13: 68QG - one number correct position, one letter wrong position
        lambda: (n1 == '6' and n2 == '8') and ((l1 == 'Q' and l2 != 'G') or (l1 != 'Q' and l2 == 'G'))
    ]
    
    return all(cond() for cond in conditions)

candidates = [['4', '8', 'S', 'Q'], ['5', '8', 'S', 'Q'], ['9', '8', 'S', 'Q']]
for pwd in candidates:
    if verify_all_conditions(pwd):
        print(pwd)