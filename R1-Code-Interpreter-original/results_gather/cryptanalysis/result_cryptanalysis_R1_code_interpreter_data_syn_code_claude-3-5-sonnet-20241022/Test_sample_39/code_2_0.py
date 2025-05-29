def validate_password(password):
    n1, n2, l1, l2 = password
    
    # For 87ID
    # 7 must be in password but in wrong position compared to 87
    if n1 == '7' and '7' == '8': return False
    if n2 == '7' and '7' == '7': return False
    if '7' not in [n1, n2]: return False
    
    # 8 cannot be in password as it's too large
    if '8' in [n1, n2]: return False
    
    # For 50FH
    # Both numbers must be larger than 5 and 0
    if int(n1) <= 5 or int(n2) <= 5: return False
    
    # Letter conditions
    # One letter from each guess must be correct but in wrong position
    # From 87ID: either I or D
    # From 50FH: either F or H
    if not ((l1 in 'ID' and l2 in 'FH') or (l2 in 'ID' and l1 in 'FH')):
        return False
    
    # One letter is too early in alphabet for each guess
    # Check if positions make sense with alphabetical order
    if l1 < l2 and (l1 not in 'FHID' or l2 not in 'FHID'):
        return False
    if l2 < l1 and (l2 not in 'FHID' or l1 not in 'FHID'):
        return False
    
    return True

valid_passwords = []
for n1 in range(6, 10):  # must be larger than 5
    for n2 in range(6, 10):
        if n1 != n2:
            for l1 in 'DFHI':  # only these letters are possible
                for l2 in 'DFHI':
                    if l1 != l2:
                        password = [str(n1), str(n2), l1, l2]
                        if validate_password(password):
                            valid_passwords.append(password)

print(valid_passwords)