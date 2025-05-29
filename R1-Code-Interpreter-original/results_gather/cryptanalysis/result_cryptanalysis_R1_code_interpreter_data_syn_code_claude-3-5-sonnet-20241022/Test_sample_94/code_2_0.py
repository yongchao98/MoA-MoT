from itertools import product
import string

def is_valid_password(password, guesses_feedback):
    for guess, feedback in guesses_feedback:
        pwd = list(password)
        gss = list(guess)
        
        # For 61VG guess
        if guess == "61VG":
            # 6 must be in first position
            if pwd[0] != '6':
                return False
            # 1 is too small compared to second number
            if pwd[1] <= '1':
                return False
            # G must be in third position
            if pwd[2] != 'G':
                return False
            # The last letter must be before V
            if pwd[3] >= 'V':
                return False
            # One letter is correct but in wrong position (G)
            # Other letter must be incorrect and too late in alphabet
            if pwd[3] >= 'V':
                return False
        
        # For guesses with "both letters are incorrect"
        if "both letters are incorrect" in feedback:
            if gss[2] in pwd[2:] or gss[3] in pwd[2:]:
                return False
            
        # For 40SV, 57UY, 58RV - both letters too late in alphabet
        if guess in ["40SV", "57UY", "58RV"]:
            # If any letter in guess is before or equal to our password letters
            if gss[2] <= pwd[2] or gss[3] <= pwd[3]:
                return False
                
        # For 01ZD - both numbers too small
        if guess == "01ZD":
            if int(gss[0]) >= int(pwd[0]) or int(gss[1]) >= int(pwd[1]):
                return False
                
        # For 87HB, 87FM, 97SH - both numbers too large
        if guess in ["87HB", "87FM", "97SH"]:
            if int(gss[0]) <= int(pwd[0]) or int(gss[1]) <= int(pwd[1]):
                return False
                
        # For 38BY
        if guess == "38BY":
            # 3 must be in second position (since we know 6 is first)
            if pwd[1] != '3':
                return False
                
        # For specific letter sequence checking
        if guess == "47QO":
            if 'Q' == pwd[3] or 'O' == pwd[3]:
                return False
        if guess == "94WK":
            if 'W' == pwd[3] or 'K' == pwd[3]:
                return False
        if guess == "29TE":
            if 'T' == pwd[3] or 'E' == pwd[3]:
                return False
        if guess == "18HS":
            if 'H' == pwd[3] or 'S' == pwd[3]:
                return False
        if guess == "92LU":
            if 'L' == pwd[3] or 'U' == pwd[3]:
                return False

    # Additional constraint: The last letter must be after M (from 68HM feedback)
    if pwd[3] <= 'M':
        return False
    
    # Additional constraint: The last letter must be before R (from 58RV feedback)
    if pwd[3] >= 'R':
        return False

    return True

# Define all guesses and their feedback (same as before)
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

valid_passwords = []
for l2 in string.ascii_uppercase:
    if l2 != 'G':  # G must be in third position
        password = '6' + '3' + 'G' + l2
        if is_valid_password(password, guesses_feedback):
            valid_passwords.append(['6', '3', 'G', l2])

print(valid_passwords)