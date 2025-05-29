def verify_password(password):
    # Convert password to string format for easier comparison
    pwd = ''.join(map(str, password))
    
    # Clue 1: 14HF - both numbers wrong, one letter correct wrong pos, one wrong
    if '1' in pwd or '4' in pwd:
        return False
    if pwd.count('H') + pwd.count('F') != 1:  # exactly one of H or F
        return False
    
    # Clue 2: 38FX - one number right wrong pos, one too large, both letters wrong
    if 'F' in pwd or 'X' in pwd:
        return False
    
    # Clue 3: 52YF - one number right wrong pos, one too large, both letters wrong
    if 'Y' in pwd:
        return False
    
    # Clue 4: 32QI - both numbers right wrong pos, letters too late
    if not (pwd[0:2] == '32' or pwd[0:2] == '23'):
        return False
    if 'Q' in pwd or 'I' in pwd:
        return False
    
    # Clue 5: 71DX - both numbers wrong, one letter right pos, one too late
    if '7' in pwd or '1' in pwd:
        return False
    if 'D' not in pwd:
        return False
    
    return True

# Test the specific password we think is correct
password = ['3', '2', 'D', 'H']
if verify_password(password):
    print(password)