def validate_password(pwd):
    def is_before(a, b):
        return ord(a) < ord(b)
    
    # Helper function to check if a number is too small
    def is_too_small(test, target):
        return int(test) < int(target)
    
    # Check 74JY
    def check_74JY():
        # Count conditions
        correct_pos_num = 0
        small_num = 0
        wrong_pos_letter = 0
        early_letter = 0
        
        # Check numbers
        if pwd[0] == '7' or pwd[1] == '4':
            correct_pos_num += 1
        if is_too_small(pwd[0], '7') or is_too_small(pwd[1], '4'):
            small_num += 1
            
        # Check letters
        if 'J' in pwd[2:] or 'Y' in pwd[2:]:
            wrong_pos_letter += 1
        if is_before('J', pwd[2]) or is_before('Y', pwd[3]):
            early_letter += 1
            
        return (correct_pos_num == 1 and small_num == 1 and 
                wrong_pos_letter == 1 and early_letter == 1)
    
    # Check 93ZN
    def check_93ZN():
        return not ('9' in pwd or '3' in pwd or 'Z' in pwd or 'N' in pwd)
    
    # Check 26MU
    def check_26MU():
        return (is_too_small('2', pwd[0]) and is_too_small('2', pwd[1]) and
                is_too_small('6', pwd[0]) and is_too_small('6', pwd[1]) and
                'M' not in pwd and 'U' not in pwd)
    
    # Check 57FS
    def check_57FS():
        # Count conditions
        wrong_pos_num = 0
        small_num = 0
        correct_pos_letter = 0
        early_letter = 0
        
        # Check numbers
        if '5' in pwd[1:] or '7' in pwd and pwd[0] != '7':
            wrong_pos_num += 1
        if is_too_small('5', pwd[0]) or is_too_small('5', pwd[1]):
            small_num += 1
            
        # Check letters
        if pwd[2] == 'F' or pwd[3] == 'S':
            correct_pos_letter += 1
        if is_before('F', pwd[2]) or is_before('S', pwd[3]):
            early_letter += 1
            
        return (wrong_pos_num == 1 and small_num == 1 and 
                correct_pos_letter == 1 and early_letter == 1)
    
    return (check_74JY() and check_93ZN() and 
            check_26MU() and check_57FS())

# Generate all possible combinations
valid_passwords = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        pwd = n1 + n2 + l1 + l2
                        if validate_password(pwd):
                            valid_passwords.append(pwd)

print(valid_passwords)