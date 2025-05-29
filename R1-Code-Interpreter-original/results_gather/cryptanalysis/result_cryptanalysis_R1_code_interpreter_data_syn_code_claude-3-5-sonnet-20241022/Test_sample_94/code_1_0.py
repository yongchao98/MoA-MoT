from itertools import product
import string

def is_valid_password(password, guesses_feedback):
    for guess, feedback in guesses_feedback:
        # Convert password and guess to lists for easier comparison
        pwd = list(password)
        gss = list(guess)
        
        # Check specific conditions for each guess
        if guess == "38BY":
            # one number is correct but in wrong position
            if not (('3' in pwd[:2] and pwd.index('3') != 0) or 
                   ('8' in pwd[:2] and pwd.index('8') != 1)):
                return False
            # one number incorrect and too large
            if not ('8' > pwd[pwd.index('3')] if '3' in pwd[:2] else '3' > pwd[1]):
                return False
        
        elif guess == "68HM":
            # one number correct and in correct position
            if pwd[0] != '6':  # Based on multiple guesses, 6 must be in first position
                return False
            
        elif guess == "61VG":
            # 6 must be in correct position (first position)
            if pwd[0] != '6':
                return False
            # 1 is too small compared to second number
            if pwd[1] <= '1':
                return False
            # G must be in third position (based on letter patterns)
            if pwd[2] != 'G':
                return False
            
        # Check number conditions
        if "both numbers are incorrect and too large" in feedback:
            for i in range(2):
                if int(gss[i]) <= int(pwd[i]):
                    return False
                    
        if "both numbers are incorrect and too small" in feedback:
            for i in range(2):
                if int(gss[i]) >= int(pwd[i]):
                    return False
                    
        # Check letter conditions
        if "both letters are incorrect and too late in the alphabet" in feedback:
            for i in range(2, 4):
                if gss[i] <= pwd[i]:
                    return False
                    
        if "both letters are incorrect" in feedback:
            for i in range(2, 4):
                if gss[i] == pwd[i]:
                    return False

    return True

# Define all guesses and their feedback
guesses_feedback = [
    ("57UY", "both numbers are incorrect; both letters are incorrect and too late in the alphabet"),
    ("40SV", "both numbers are incorrect; both letters are incorrect and too late in the alphabet"),
    ("01ZD", "both numbers are incorrect and too small; both letters are incorrect"),
    ("58RV", "both numbers are incorrect; both letters are incorrect and too late in the alphabet"),
    ("90FV", "both numbers are incorrect; both letters are incorrect"),
    ("38BY", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect"),
    ("47QO", "both numbers are incorrect; both letters are incorrect"),
    ("94WK", "both numbers are incorrect; both letters are incorrect"),
    ("68HM", "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect"),
    ("87HB", "both numbers are incorrect and too large; both letters are incorrect"),
    ("97SH", "both numbers are incorrect and too large; both letters are incorrect"),
    ("29TE", "both numbers are incorrect; both letters are incorrect"),
    ("18HS", "both numbers are incorrect; both letters are incorrect"),
    ("92LU", "both numbers are incorrect; both letters are incorrect"),
    ("87FM", "both numbers are incorrect and too large; both letters are incorrect"),
    ("61VG", "one number is correct and in the correct position; one number is incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and too late in the alphabet")
]

# Generate all possible combinations
numbers = list(map(str, range(10)))
letters = list(string.ascii_uppercase)
valid_passwords = []

# Based on the feedback, we know:
# 1. First number must be 6
# 2. Second number must be 3 (from multiple conditions)
# 3. Third position must be G (from 61VG feedback and other constraints)
# 4. Last letter must be before V in alphabet

for l2 in letters:
    if l2 != 'G' and l2 < 'V':  # Must be before V
        password = '6' + '3' + 'G' + l2
        if is_valid_password(password, guesses_feedback):
            valid_passwords.append(['6', '3', 'G', l2])

print(valid_passwords)