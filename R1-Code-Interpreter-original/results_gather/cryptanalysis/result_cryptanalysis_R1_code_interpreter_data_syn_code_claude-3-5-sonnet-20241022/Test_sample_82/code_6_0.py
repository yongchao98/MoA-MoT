def check_conditions(n1, n2, l1, l2):
    # Immediate rejection if any digits or letters repeat
    if n1 == n2 or l1 == l2:
        return False
        
    # Convert to list for easier handling
    password = [n1, n2, l1, l2]
    
    # Condition 1: 87CF
    if int(n1) >= 8 or int(n2) >= 7:  # numbers must be smaller
        return False
    if ord(l1) <= ord('C') or ord(l2) <= ord('F'):  # letters must be later
        return False
    
    # Condition 2: 71SM
    # Exactly one of 7 or 1 must be correct but in wrong position
    if '7' in [n1, n2]:  # 7 should not be in answer as it's too large
        return False
    if [n1, n2].count('1') != 1:  # exactly one 1 should be present
        return False
    
    # One of S or M must be correct but in wrong position
    if ('S' in [l1, l2] and 'M' in [l1, l2]):  # can't have both
        return False
    if not ('S' in [l1, l2] or 'M' in [l1, l2]):  # must have one
        return False
    
    # Condition 3: 23AY
    if '2' in [n1, n2] or '3' in [n1, n2]:
        return False
    if 'A' in [l1, l2] or 'Y' in [l1, l2]:
        return False
    
    # Condition 4: 53PD
    if '5' in [n1, n2] or '3' in [n1, n2]:
        return False
    # One of P or D must be correct but in wrong position
    if ('P' in [l1, l2] and 'D' in [l1, l2]):  # can't have both
        return False
    if not ('P' in [l1, l2] or 'D' in [l1, l2]):  # must have one
        return False
    
    # Condition 5: 01UF
    # Check for 0 and 1
    zero_one_count = [n1, n2].count('0') + [n1, n2].count('1')
    if zero_one_count != 1:  # exactly one of 0 or 1
        return False
    if '0' in [n1, n2]:  # if 0 present, it must be too small
        if int(n1) == 0 or int(n2) == 0:
            return False
    if 'U' in [l1, l2] or 'F' in [l1, l2]:
        return False
    
    # Condition 6: 43US
    if '4' in [n1, n2] or '3' in [n1, n2]:
        return False
    if not (ord(l1) < ord('U') and ord(l2) < ord('S')):  # letters must be earlier
        return False
    
    return True

def find_password():
    valid_numbers = [str(i) for i in range(10)]
    valid_letters = [chr(i) for i in range(ord('G'), ord('R'))]  # Reasonable range
    
    for n1 in valid_numbers:
        for n2 in valid_numbers:
            if n1 == n2:  # skip if numbers are same
                continue
            for l1 in valid_letters:
                for l2 in valid_letters:
                    if l1 == l2:  # skip if letters are same
                        continue
                    if check_conditions(n1, n2, l1, l2):
                        return [n1, n2, l1, l2]
    return None

result = find_password()
print(result)