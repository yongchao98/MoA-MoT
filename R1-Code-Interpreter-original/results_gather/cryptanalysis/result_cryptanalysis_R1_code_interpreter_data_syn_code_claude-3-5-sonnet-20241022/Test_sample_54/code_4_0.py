def check_guess(candidate):
    n1, n2, l1, l2 = candidate
    
    # Condition 1: 04XW - both numbers wrong, letters too late
    if n1 == '0' or n2 == '4' or l1 >= 'X' or l2 >= 'W':
        return False
    
    # Condition 2: 06IH - numbers wrong, H must be present but in WRONG position
    if n1 == '0' or n2 == '6':
        return False
    if 'H' not in [l1, l2]:  # H must be one of the letters
        return False
    if l2 == 'H':  # H must be in wrong position compared to 06IH
        return False
    
    # Condition 3: 54FA - all incorrect
    if n1 == '5' or n2 == '4' or l1 == 'F' or l2 == 'A':
        return False
    
    # Condition 4: 07ES - numbers wrong, E must be present but in WRONG position
    if n1 == '0' or n2 == '7':
        return False
    if 'E' not in [l1, l2]:  # E must be one of the letters
        return False
    if l2 == 'E':  # E must be in wrong position compared to 07ES
        return False
    
    # Condition 5: 08HI - one number (8) correct in pos
    if n2 != '8':  # 8 must be in second position
        return False
    
    # Condition 6: 57OR - numbers wrong, letters too late
    if n1 == '5' or n2 == '7' or l1 >= 'O' or l2 >= 'O':
        return False
    
    # Condition 7: 29FB - one number (2) correct in pos
    if n1 != '2':  # 2 must be in first position
        return False
    if l1 == 'F' or l2 == 'B':  # Letters can't be F or B
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        if check_guess([n1, n2, l1, l2]):
                            valid_solutions.append([n1, n2, l1, l2])

print(valid_solutions)