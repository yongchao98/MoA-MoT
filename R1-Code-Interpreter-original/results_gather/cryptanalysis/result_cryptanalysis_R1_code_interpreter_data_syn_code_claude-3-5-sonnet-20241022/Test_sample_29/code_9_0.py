def validate_guess(pwd):
    def is_before(a, b):
        return ord(a) < ord(b)
    
    # Basic validation
    if len(pwd) != 4 or pwd[0] == pwd[1] or pwd[2] == pwd[3]:
        return False
        
    # 74JY conditions
    c1_num_correct = (pwd[1] == '4')  # 4 must be in second position
    c1_num_small = int(pwd[0]) < 7  # first number must be smaller than 7
    c1_letter_wrong = ('Y' in pwd[2:])  # Y must be in password but wrong position
    c1_letter_early = is_before('J', pwd[2] if pwd[2] != 'J' else pwd[3])
    
    if not (c1_num_correct and c1_num_small and c1_letter_wrong and c1_letter_early):
        return False
    
    # 93ZN conditions - none of these can be in password
    if '9' in pwd or '3' in pwd or 'Z' in pwd or 'N' in pwd:
        return False
    
    # 26MU conditions
    if ('2' in pwd or '6' in pwd or 'M' in pwd or 'U' in pwd):
        return False
    if not (int(pwd[0]) > 2 and int(pwd[1]) > 2):  # both numbers must be > 2
        return False
    
    # 57FS conditions
    c4_num_wrong = ('5' in pwd and pwd[1] != '5')  # 5 must be in password but not in position 2
    c4_num_small = int(pwd[0]) < 5 or int(pwd[1]) < 7  # one number must be smaller
    c4_letter_correct = (pwd[3] == 'S')  # S must be in last position
    c4_letter_early = is_before('F', pwd[2])  # F must be before the other letter
    
    return (c4_num_wrong and c4_num_small and c4_letter_correct and c4_letter_early)

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
                        if validate_guess(pwd):
                            valid_passwords.append(pwd)

print(valid_passwords)