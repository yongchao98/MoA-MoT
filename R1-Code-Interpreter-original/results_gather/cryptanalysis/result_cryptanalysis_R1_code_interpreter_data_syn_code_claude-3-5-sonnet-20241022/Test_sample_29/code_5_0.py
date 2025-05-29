def check_conditions(pwd):
    # Basic format check
    if len(pwd) != 4 or not (pwd[0].isdigit() and pwd[1].isdigit() and 
                            pwd[2].isalpha() and pwd[3].isalpha()):
        return False
    
    # Check for no repeats
    if pwd[0] == pwd[1] or pwd[2] == pwd[3]:
        return False
    
    # Condition 1: 74JY
    c1_num_correct = (pwd[0] == '4')  # 4 must be in first position
    c1_num_small = int(pwd[1]) > 4    # other number must be larger than 4
    c1_letter_wrong = ('Y' in pwd[2:])  # Y must be in password but wrong position
    c1_letter_early = (ord('J') < ord(pwd[2]) and ord('J') < ord(pwd[3]))  # J too early
    
    if not (c1_num_correct and c1_num_small and c1_letter_wrong and c1_letter_early):
        return False
    
    # Condition 2: 93ZN
    if '9' in pwd or '3' in pwd or 'Z' in pwd or 'N' in pwd:
        return False
    
    # Condition 3: 26MU
    if ('2' in pwd or '6' in pwd or 'M' in pwd or 'U' in pwd or 
        int(pwd[0]) < 2 or int(pwd[1]) < 6):
        return False
    
    # Condition 4: 57FS
    c4_num_wrong = ('5' in pwd and pwd[0] != '5')  # 5 must be in password but wrong position
    c4_num_small = True  # will be checked separately
    c4_letter_correct = ('S' == pwd[3])  # S must be in correct position (last)
    c4_letter_early = (ord('F') < ord(pwd[2]))  # F too early
    
    if not (c4_num_wrong and c4_letter_correct and c4_letter_early):
        return False
    
    return True

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
                        if check_conditions(pwd):
                            valid_passwords.append(pwd)

print(valid_passwords)