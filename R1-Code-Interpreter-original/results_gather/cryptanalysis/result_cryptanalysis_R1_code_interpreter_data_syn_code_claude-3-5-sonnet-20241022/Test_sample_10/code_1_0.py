def is_valid_password(password):
    # Convert password tuple to list for easier handling
    pwd = list(password)
    
    # Guess 1: 25KS - all wrong
    if '2' in pwd or '5' in pwd or 'K' in pwd or 'S' in pwd:
        return False
    
    # Guess 2: 26BP
    if '2' in pwd or '6' in pwd:
        return False
    # One letter (B/P) correct but wrong position, other too late
    if not ((('B' in pwd and pwd.index('B') != 3) and ('P' > pwd[2] and 'P' > pwd[3])) or
            (('P' in pwd and pwd.index('P') != 3) and ('B' > pwd[2] and 'B' > pwd[3]))):
        return False
    
    # Guess 3: 39ON
    # One number (3/9) correct but wrong position, other too small
    if not ((('3' in pwd and pwd.index('3') != 0) and int('9') < int(pwd[0])) or
            (('9' in pwd and pwd.index('9') != 0) and int('3') < int(pwd[0]))):
        return False
    # One letter (O/N) correct but wrong position, other too late
    if not ((('O' in pwd and pwd.index('O') != 3) and ('N' > pwd[2] and 'N' > pwd[3])) or
            (('N' in pwd and pwd.index('N') != 3) and ('O' > pwd[2] and 'O' > pwd[3]))):
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10) if i not in [2, 5, 6]]
letters = [c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' if c not in ['K', 'S']]

valid_passwords = []
from itertools import permutations

for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        pwd = (n1, n2, l1, l2)
        if is_valid_password(pwd):
            valid_passwords.append(pwd)

print(valid_passwords)