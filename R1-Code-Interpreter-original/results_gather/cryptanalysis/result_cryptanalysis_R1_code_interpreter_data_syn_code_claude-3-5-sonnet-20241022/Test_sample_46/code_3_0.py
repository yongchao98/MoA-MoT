def verify_guess(password):
    # Guess 3 (09OQ)
    if password[0] != '0' or password[1] != '9':  # Numbers must be 09
        return False
        
    # From guess 3: One of O/Q must be in password but in wrong position
    if 'O' in password[2:]:
        if password.index('O') == 2:  # O can't be in first letter position
            return False
    elif 'Q' in password[2:]:
        if password.index('Q') == 3:  # Q can't be in second letter position
            return False
    else:
        return False
        
    # From guess 4: R must be in password but in wrong position
    if 'R' not in password[2:]:
        return False
    if password.index('R') == 3:  # R can't be in second letter position
        return False
        
    # Letters must be after M in alphabet
    if any(c <= 'M' for c in password[2:]):
        return False
        
    # Letters can't repeat
    if password[2] == password[3]:
        return False
        
    return True

# Generate all possible combinations
valid_passwords = []
numbers = '09'
import string
letters = [l for l in string.ascii_uppercase if l > 'M']

for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            password = ['0', '9', l1, l2]
            if verify_guess(password):
                valid_passwords.append(password)

print(valid_passwords)