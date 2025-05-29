def is_valid_password(pwd):
    # Helper function to check if a letter comes before another in alphabet
    def comes_before(a, b):
        return ord(a) < ord(b)
    
    # Condition 1: 74JY
    def check_74JY():
        correct_pos_num = (pwd[0] == '7' or pwd[1] == '4')
        other_num_small = (pwd[0] < '7' or pwd[1] < '4')
        letter_wrong_pos = ('J' in pwd[2:] or 'Y' in pwd[2:])
        other_letter_early = (comes_before('J', pwd[2]) or comes_before('Y', pwd[3]))
        return correct_pos_num and other_num_small and letter_wrong_pos and other_letter_early
    
    # Condition 2: 93ZN
    def check_93ZN():
        return not ('9' in pwd or '3' in pwd or 'Z' in pwd or 'N' in pwd)
    
    # Condition 3: 26MU
    def check_26MU():
        nums_too_small = int(pwd[0]) > 2 and int(pwd[0]) > 6 and int(pwd[1]) > 2 and int(pwd[1]) > 6
        letters_not_present = 'M' not in pwd and 'U' not in pwd
        return nums_too_small and letters_not_present
    
    # Condition 4: 57FS
    def check_57FS():
        num_wrong_pos = ('5' in pwd or '7' in pwd) and pwd[0] != '5' and pwd[1] != '7'
        other_num_small = int(pwd[0]) > 5 or int(pwd[1]) > 7
        letter_correct_pos = (pwd[2] == 'F' or pwd[3] == 'S')
        other_letter_early = comes_before('F', pwd[2]) or comes_before('S', pwd[3])
        return num_wrong_pos and other_num_small and letter_correct_pos and other_letter_early
    
    # Additional checks
    def basic_checks():
        # Numbers must be different
        if pwd[0] == pwd[1]:
            return False
        # Letters must be different
        if pwd[2] == pwd[3]:
            return False
        return True
    
    return (check_74JY() and check_93ZN() and check_26MU() and 
            check_57FS() and basic_checks())

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
                        if is_valid_password(pwd):
                            valid_passwords.append(pwd)

print(valid_passwords)