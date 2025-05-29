def is_valid_password(password):
    n1, n2, l1, l2 = password
    
    # Guess 1: 87CF
    if not (int('8') > int(n1) and int('7') > int(n2)):
        return False
    if not (ord('C') < ord(l1) and ord('F') < ord(l2)):
        return False
    
    # Guess 2: 71SM
    num_correct = 0
    if '7' in [n1, n2]: num_correct += 1
    if '1' in [n1, n2]: num_correct += 1
    if num_correct != 1:  # exactly one number correct but wrong position
        return False
    if '7' in [n1, n2] and int('7') <= max(int(n1), int(n2)):  # 7 must be too large
        return False
    
    # Important: exactly one of S or M must be in the password but in wrong position
    if ('S' == l1 or 'S' == l2) and ('M' == l1 or 'M' == l2):
        return False
    if not ('S' in [l1, l2] or 'M' in [l1, l2]):
        return False
    
    # Guess 3: 23AY
    if '2' in [n1, n2] or '3' in [n1, n2]:
        return False
    if 'A' in [l1, l2] or 'Y' in [l1, l2]:
        return False
    
    # Guess 4: 53PD
    if '5' in [n1, n2] or '3' in [n1, n2]:
        return False
    if not ('P' in [l1, l2] or 'D' in [l1, l2]):
        return False
    
    # Guess 5: 01UF
    if [n1, n2].count('0') + [n1, n2].count('1') != 1:
        return False
    if '0' in [n1, n2] and int('0') >= min(int(n1), int(n2)):  # 0 must be too small
        return False
    if 'U' in [l1, l2] or 'F' in [l1, l2]:
        return False
    
    # Guess 6: 43US
    if '4' in [n1, n2] or '3' in [n1, n2]:
        return False
    if not (ord('U') > ord(l1) and ord('S') > ord(l2)):
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