from itertools import product

def check_alphabet_order(password):
    # From guess 6: if P is not in position, L must be the one in wrong position
    # and the other letter must be later in alphabet than L
    if 'P' in password[2:]:
        other_letter = password[2] if password[3] == 'P' else password[3]
        if other_letter <= 'L':  # other letter must be later than L
            return False
    if 'L' in password[2:]:
        other_letter = password[2] if password[3] == 'L' else password[3]
        if other_letter <= 'L':  # other letter must be later than L
            return False
    return True

def is_valid_password(password):
    # Convert password to list of strings
    password = [str(x) for x in password]
    
    # Basic conditions from previous code
    if password[0] == password[1]: return False  # numbers can't repeat
    if password[2] == password[3]: return False  # letters can't repeat
    
    # Number conditions
    if password[0] != '7': return False  # 7 must be in first position
    if password[1] != '3': return False  # 3 must be in second position
    
    # Letter conditions
    if not ('Q' in password[2:] or 'Y' in password[2:]): return False
    if not ('P' in password[2:] or 'L' in password[2:]): return False
    
    # Check alphabet order condition
    if not check_alphabet_order(password):
        return False
    
    return True

# Generate all possible combinations
numbers = [7, 3]  # We know these are fixed
letters = ['L', 'P', 'Q', 'Y']  # These are our possible letters
valid_passwords = []

for combo in product(numbers, numbers, letters, letters):
    if is_valid_password(combo):
        valid_passwords.append(combo)

print(valid_passwords)