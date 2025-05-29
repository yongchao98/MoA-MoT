def verify_guess(pwd):
    # Helper function to check if a letter comes before another in alphabet
    def is_before(a, b):
        return ord(a) < ord(b)
    
    # Guess 1: 74JY
    def check_74JY():
        # Count correct numbers and positions
        nums_correct_pos = (pwd[0] == '7' or pwd[1] == '4')
        # Count numbers that are too small
        nums_too_small = (int(pwd[0]) < 7 or int(pwd[1]) < 4)
        # Count letters in wrong position
        letters_wrong_pos = ('J' in pwd[2:] or 'Y' in pwd[2:])
        # Count letters that are too early
        letters_early = (is_before('J', pwd[2]) or is_before('Y', pwd[3]))
        
        return (nums_correct_pos and nums_too_small and 
                letters_wrong_pos and letters_early)
    
    # Guess 2: 93ZN
    def check_93ZN():
        return not ('9' in pwd or '3' in pwd or 
                   'Z' in pwd or 'N' in pwd)
    
    # Guess 3: 26MU
    def check_26MU():
        nums_too_small = (int(pwd[0]) > 2 and int(pwd[1]) > 2 and 
                         int(pwd[0]) > 6 and int(pwd[1]) > 6)
        no_letters = ('M' not in pwd and 'U' not in pwd)
        return nums_too_small and no_letters
    
    # Guess 4: 57FS
    def check_57FS():
        # One number correct but wrong position
        nums_wrong_pos = (('5' in pwd and pwd.index('5') != 0) or 
                         ('7' in pwd and pwd.index('7') != 1))
        # One number too small
        nums_too_small = (int(pwd[0]) < 5 or int(pwd[1]) < 7)
        # One letter correct and in position
        letters_correct_pos = (pwd[2] == 'F' or pwd[3] == 'S')
        # One letter too early
        letters_early = (is_before('F', pwd[2]) or is_before('S', pwd[3]))
        
        return (nums_wrong_pos and nums_too_small and 
                letters_correct_pos and letters_early)
    
    return (check_74JY() and check_93ZN() and 
            check_26MU() and check_57FS())

# Generate all possible combinations
valid_passwords = []
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        pwd = n1 + n2 + l1 + l2
                        if verify_guess(pwd):
                            valid_passwords.append(pwd)

print(valid_passwords)