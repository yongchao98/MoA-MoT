def check_conditions(password):
    # Convert password to list of strings for easier comparison
    pwd = [str(x) for x in password]
    
    # Condition 1: 14LZ
    # Both numbers (1,4) correct but wrong positions
    if not (set(pwd[:2]) == set(['1', '4']) and 
            pwd[0] != '1' and pwd[1] != '4'):
        return False
    
    # Condition 2: 28EP
    # Both numbers incorrect, one letter correct in position, one letter too early
    if ('2' in pwd[:2] or '8' in pwd[:2]):
        return False
    # Check if P is in position 4 (index 3) - this must be true because E is too early
    p_correct = (pwd[3] == 'P')
    # The other letter in pwd must be later than E
    other_letter_pos = 2
    if p_correct and ord(pwd[other_letter_pos]) <= ord('E'):
        return False
    
    # Condition 3: 14RW
    # Both numbers correct but wrong positions (already checked in condition 1)
    if 'R' in pwd[2:] or 'W' in pwd[2:]:
        return False
    
    # Condition 4: 70XH
    # Both numbers incorrect, one letter correct in position, one letter too early
    if ('7' in pwd[:2] or '0' in pwd[:2]):
        return False
    # Since P is in position 4 from condition 2, H cannot be in position 4
    if pwd[3] == 'H':
        return False
    # The letter in position 3 must be later than X
    if ord(pwd[2]) <= ord('X'):
        return False
    
    # Additional checks:
    # Letters can't be L,Z,R,W,E,X
    if any(x in pwd[2:] for x in ['L','Z','R','W','E','X']):
        return False
    
    return True

# Generate all possible combinations
numbers = ['4', '1']  # We know these are the correct numbers
letters = [chr(i) for i in range(ord('A'), ord('Z')+1)]
letters = [l for l in letters if l not in 'LZRWEXH']  # Remove known incorrect letters

valid_passwords = []
for n1 in ['1', '4']:
    n2 = '4' if n1 == '1' else '1'
    for l1 in letters:
        if ord(l1) > ord('X'):  # First letter must be after X
            for l2 in ['P']:  # Second letter must be P
                if l1 != l2:  # Letters can't repeat
                    password = [n1, n2, l1, l2]
                    if check_conditions(password):
                        valid_passwords.append(password)

print(valid_passwords)