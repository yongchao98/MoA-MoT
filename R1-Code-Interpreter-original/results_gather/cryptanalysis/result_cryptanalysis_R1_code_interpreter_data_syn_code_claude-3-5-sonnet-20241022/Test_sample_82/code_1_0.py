def is_valid_password(password):
    # Helper function to check if a number is too large/small compared to target
    def compare_numbers(guess, target):
        if guess == target:
            return "correct"
        return "too large" if int(guess) > int(target) else "too small"
    
    # Helper function to check if a letter is too early/late compared to target
    def compare_letters(guess, target):
        if guess == target:
            return "correct"
        return "too late" if guess > target else "too early"
    
    # Check Guess 1: 87CF
    if not (int('8') > int(password[0]) and int('7') > int(password[1])):  # both numbers too large
        return False
    if not (ord('C') < ord(password[2]) and ord('F') < ord(password[3])):  # both letters too early
        return False
    
    # Check Guess 2: 71SM
    correct_num = False
    if password[0] == '1' or password[1] == '7':  # one number correct but wrong position
        correct_num = True
    if not correct_num:
        return False
    if not (max(password[0], password[1]) < '7'):  # one number too large
        return False
    correct_letter = False
    if 'S' in password[2:] or 'M' in password[2:]:  # one letter correct but wrong position
        correct_letter = True
    if not correct_letter:
        return False
    
    # Check Guess 3: 23AY
    if '2' in password[:2] or '3' in password[:2]:  # both numbers incorrect
        return False
    if 'A' in password[2:] or 'Y' in password[2:]:  # both letters incorrect
        return False
    
    # Check Guess 4: 53PD
    if '5' in password[:2] or '3' in password[:2]:  # both numbers incorrect
        return False
    if not ('P' in password[2:] or 'D' in password[2:]):  # one letter correct but wrong position
        return False
    
    # Check Guess 5: 01UF
    correct_num = False
    if '0' in password[:2] or '1' in password[:2]:  # one number correct but wrong position
        correct_num = True
    if not correct_num:
        return False
    if 'U' in password[2:] or 'F' in password[2:]:  # both letters incorrect
        return False
    
    # Check Guess 6: 43US
    if '4' in password[:2] or '3' in password[:2]:  # both numbers incorrect
        return False
    if not (ord('U') > ord(password[2]) and ord('S') > ord(password[3])):  # both letters too late
        return False
    
    return True

# Try all possible combinations
def find_password():
    for n1 in range(10):
        for n2 in range(10):
            if n1 != n2:
                for l1 in [chr(i) for i in range(65, 91)]:  # A-Z
                    for l2 in [chr(i) for i in range(65, 91)]:  # A-Z
                        if l1 != l2:
                            password = [str(n1), str(n2), l1, l2]
                            if is_valid_password(password):
                                return password
    return None

result = find_password()
print(result)