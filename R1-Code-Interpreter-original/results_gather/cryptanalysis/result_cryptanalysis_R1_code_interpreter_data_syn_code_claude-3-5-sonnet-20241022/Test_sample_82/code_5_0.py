def check_conditions(n1, n2, l1, l2):
    # Convert to list for easier handling
    password = [n1, n2, l1, l2]
    
    # Condition 1: 87CF
    if int(n1) >= 8 or int(n2) >= 7:  # numbers must be smaller
        return False
    if ord(l1) <= ord('C') or ord(l2) <= ord('F'):  # letters must be later
        return False
    
    # Condition 2: 71SM
    has_correct_num = False
    if '1' in password[:2] and '7' not in password[:2]:
        has_correct_num = True
    if not has_correct_num:
        return False
    
    has_correct_letter = False
    if ('S' in password[2:] or 'M' in password[2:]):
        has_correct_letter = True
    if not has_correct_letter:
        return False
    
    # Condition 3: 23AY
    if '2' in password[:2] or '3' in password[:2]:
        return False
    if 'A' in password[2:] or 'Y' in password[2:]:
        return False
    
    # Condition 4: 53PD
    if '5' in password[:2] or '3' in password[:2]:
        return False
    has_correct_letter = False
    if 'P' in password[2:] or 'D' in password[2:]:
        has_correct_letter = True
    if not has_correct_letter:
        return False
    
    # Condition 5: 01UF
    has_correct_num = False
    if '1' in password[:2] or '0' in password[:2]:
        has_correct_num = True
    if not has_correct_num:
        return False
    if min(int(n1), int(n2)) <= 0:  # one number too small
        return False
    if 'U' in password[2:] or 'F' in password[2:]:
        return False
    
    # Condition 6: 43US
    if '4' in password[:2] or '3' in password[:2]:
        return False
    if ord(l1) >= ord('U') or ord(l2) >= ord('S'):  # letters must be earlier
        return False
    
    return True

def find_password():
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:  # skip if numbers are same
                continue
            n1, n2 = str(n1), str(n2)
            for l1 in [chr(i) for i in range(ord('G'), ord('S'))]:  # Reasonable letter range
                for l2 in [chr(i) for i in range(ord('G'), ord('S'))]:
                    if l1 == l2:  # skip if letters are same
                        continue
                    if check_conditions(n1, n2, l1, l2):
                        return [n1, n2, l1, l2]
    return None

result = find_password()
print(result)