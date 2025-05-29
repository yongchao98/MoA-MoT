def test_password(password):
    # Convert password list to string for easier comparison
    pwd = ''.join(password)
    
    # Guess 3 (09OQ): Both numbers correct, one letter correct but wrong pos
    if pwd[0:2] != '09':
        return False
    if not (('O' in pwd[2:] and pwd.index('O') != pwd.index('O')) or 
            ('Q' in pwd[2:] and pwd.index('Q') != pwd.index('Q'))):
        return False
    
    # Guess 4 (28BR): One letter correct but wrong pos, one too early
    if not ('R' in pwd[2:] or 'B' in pwd[2:]):
        return False
    
    # Letters must be after B,H,I,M in alphabet
    if any(l <= 'M' for l in pwd[2:]):
        return False
    
    # Letters can't repeat
    if pwd[2] == pwd[3]:
        return False
    
    return True

# We know numbers are 09
numbers = '09'
valid_passwords = []

# Try all possible letter combinations
import string
letters = [l for l in string.ascii_uppercase if l > 'M']  # Only letters after M

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # Letters can't repeat
            password = ['0', '9', l1, l2]
            if test_password(password):
                valid_passwords.append(password)

print(valid_passwords)