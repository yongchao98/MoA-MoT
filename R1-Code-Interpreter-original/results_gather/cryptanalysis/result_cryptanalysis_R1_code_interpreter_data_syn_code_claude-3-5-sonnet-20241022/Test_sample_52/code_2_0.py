def check_guess(password):
    # Convert password to list of characters for easier handling
    pwd = list(password)
    
    # Guess 4: 06GE - This is the most informative guess
    # We know one number (0 or 6) is correct and in position
    # We know one letter (G or E) is correct and in position
    if (pwd[0] != '0' and pwd[1] != '6') or (pwd[2] != 'G' and pwd[3] != 'E'):
        return False
    if pwd.count('0') + pwd.count('6') != 1 or pwd.count('G') + pwd.count('E') != 1:
        return False

    # Guess 8: 15IG - G must be in password but in wrong position
    # If G is in password, it must be in wrong position when guessed in position 3
    if 'G' in pwd and pwd[2] == 'G':
        return False

    # Guess 10: 93PA - one number (9 or 3) must be in password but in wrong position
    if not (('9' in pwd[:2] and pwd[0] != '9' and pwd[1] != '9') or 
            ('3' in pwd[:2] and pwd[0] != '3' and pwd[1] != '3')):
        return False

    # Check that numbers that must not appear don't appear
    forbidden_nums = {'1', '2', '5', '7'}
    if any(n in forbidden_nums for n in pwd[:2]):
        return False

    # Check that letters that must not appear don't appear
    forbidden_letters = {'N', 'T', 'O', 'D', 'B', 'I', 'F', 'Y', 'U', 'P', 'A'}
    if any(l in forbidden_letters for l in pwd[2:]):
        return False

    # Additional alphabet position checks
    # Letters N, T, U must be too late compared to actual letters
    if not (ord(pwd[2]) < ord('N') and ord(pwd[2]) < ord('T') and 
            ord(pwd[3]) < ord('N') and ord(pwd[3]) < ord('T') and
            ord(pwd[2]) < ord('U') and ord(pwd[3]) < ord('U')):
        return False

    return True

# Generate and test all possible combinations
valid_passwords = []
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2:
                        password = f"{n1}{n2}{l1}{l2}"
                        if check_guess(password):
                            valid_passwords.append([str(n1), str(n2), l1, l2])

print(valid_passwords)