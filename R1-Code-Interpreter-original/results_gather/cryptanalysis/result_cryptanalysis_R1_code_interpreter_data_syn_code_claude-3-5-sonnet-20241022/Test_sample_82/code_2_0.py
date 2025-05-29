def is_valid_password(password):
    n1, n2, l1, l2 = password
    
    # Guess 1: 87CF
    if not (int('8') > int(n1) and int('7') > int(n2)):  # both numbers too large
        return False
    if not (ord('C') < ord(l1) and ord('F') < ord(l2)):  # both letters too early
        return False
    
    # Guess 2: 71SM
    if not ((n2 == '7' or n1 == '1') and not (n1 == '7' or n2 == '1')):  # exactly one number correct but wrong position
        return False
    if not (ord('S') < ord(l2) or ord('M') < ord(l1)):  # one letter correct but wrong position, one too late
        return False
    
    # Guess 3: 23AY
    if '2' in [n1, n2] or '3' in [n1, n2]:  # both numbers incorrect
        return False
    if 'A' in [l1, l2] or 'Y' in [l1, l2]:  # both letters incorrect
        return False
    
    # Guess 4: 53PD
    if '5' in [n1, n2] or '3' in [n1, n2]:  # both numbers incorrect
        return False
    if not ('P' in [l1, l2] or 'D' in [l1, l2]):  # one letter correct but wrong position
        return False
    
    # Guess 5: 01UF
    zero_one_count = [n1, n2].count('0') + [n1, n2].count('1')
    if zero_one_count != 1:  # exactly one number correct but wrong position
        return False
    if int(n1) < 1 or int(n2) < 1:  # one number too small
        return False
    if 'U' in [l1, l2] or 'F' in [l1, l2]:  # both letters incorrect
        return False
    
    # Guess 6: 43US
    if '4' in [n1, n2] or '3' in [n1, n2]:  # both numbers incorrect
        return False
    if not (ord('U') > ord(l1) and ord('S') > ord(l2)):  # both letters too late
        return False
    
    return True

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