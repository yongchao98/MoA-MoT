def check_conditions(pwd):
    # Convert to list for easier handling
    pwd = list(pwd)
    
    # Condition 1: 25KS - all wrong
    if '2' in pwd or '5' in pwd or 'K' in pwd or 'S' in pwd:
        return False
        
    # Condition 2: 26BP
    if '2' in pwd or '6' in pwd:
        return False
    # For 26BP: one letter correct but wrong position
    bp_condition = False
    if 'B' in pwd[2:] or 'P' in pwd[2:]:
        if 'B' in pwd[2:]:
            bp_condition = True
        elif 'P' in pwd[2:]:
            bp_condition = True
    if not bp_condition:
        return False
        
    # Condition 3: 39ON
    # One number correct but wrong position
    num_condition = False
    if '3' in pwd[1:] or '9' in pwd[1:]:
        if '3' in pwd[1:]:
            # Other number (9) should be too small
            if int(pwd[0]) > 9:
                num_condition = True
        elif '9' in pwd[1:]:
            # Other number (3) should be too small
            if int(pwd[0]) > 3:
                num_condition = True
    if not num_condition:
        return False
        
    # One letter correct but wrong position
    on_condition = False
    if ('O' in pwd[2:] and pwd.index('O') != 3) or ('N' in pwd[2:] and pwd.index('N') != 3):
        on_condition = True
    if not on_condition:
        return False
    
    return True

# Generate possible combinations
numbers = [str(i) for i in range(10) if i not in [2, 5, 6]]
letters = [c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' if c not in ['K', 'S']]

from itertools import permutations

valid_passwords = []
for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(['B', 'P', 'O', 'N'], 2):  # Focus on relevant letters first
        pwd = (n1, n2, l1, l2)
        if check_conditions(pwd):
            valid_passwords.append(pwd)

print(valid_passwords)