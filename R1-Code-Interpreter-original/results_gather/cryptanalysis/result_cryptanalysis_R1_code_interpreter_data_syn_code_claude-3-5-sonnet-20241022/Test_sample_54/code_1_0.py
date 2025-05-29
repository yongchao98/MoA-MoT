def check_conditions(guess):
    # Convert guess to list of characters
    g = list(guess)
    
    # Check all conditions
    conditions = [
        # Condition 1: 04XW - all incorrect, letters too late
        lambda g: g[0] != '0' and g[1] != '4' and g[2] < 'X' and g[3] < 'W',
        
        # Condition 2: 06IH - numbers wrong, one letter correct wrong pos, one too late
        lambda g: g[0] != '0' and g[1] != '6' and 
                 ((g[2] == 'H' and g[3] != 'I') or (g[2] != 'I' and g[3] == 'H')) and
                 (g[2] < 'H' or g[3] < 'H'),
        
        # Condition 3: 54FA - all incorrect
        lambda g: g[0] != '5' and g[1] != '4' and g[2] != 'F' and g[3] != 'A',
        
        # Condition 4: 07ES - numbers wrong, one letter correct wrong pos, one too late
        lambda g: g[0] != '0' and g[1] != '7' and
                 ((g[2] == 'S' and g[3] != 'E') or (g[2] != 'E' and g[3] == 'S')) and
                 (g[2] < 'S' or g[3] < 'S'),
        
        # Condition 5: 08HI - one number correct in pos, one too small, one letter correct in pos, one too late
        lambda g: ((g[0] == '0' and g[1] != '8') or (g[0] != '0' and g[1] == '8')) and
                 ((g[0] < '8' and g[0] != '0') or (g[1] < '8' and g[1] != '8')) and
                 ((g[2] == 'H' and g[3] != 'I') or (g[2] != 'H' and g[3] == 'I')) and
                 (g[2] < 'I' or g[3] < 'I'),
        
        # Condition 6: 57OR - numbers wrong, letters too late
        lambda g: g[0] != '5' and g[1] != '7' and g[2] < 'O' and g[3] < 'R',
        
        # Condition 7: 29FB - one number correct in pos, one too large, letters wrong
        lambda g: ((g[0] == '2' and g[1] != '9') or (g[0] != '2' and g[1] == '9')) and
                 ((g[0] > '9' and g[0] != '2') or (g[1] > '9' and g[1] != '9')) and
                 g[2] != 'F' and g[3] != 'B'
    ]
    
    return all(condition(g) for condition in conditions)

# Generate and test all possible combinations
valid_solutions = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        guess = n1 + n2 + l1 + l2
                        if check_conditions(guess):
                            valid_solutions.append([n1, n2, l1, l2])

print(valid_solutions)