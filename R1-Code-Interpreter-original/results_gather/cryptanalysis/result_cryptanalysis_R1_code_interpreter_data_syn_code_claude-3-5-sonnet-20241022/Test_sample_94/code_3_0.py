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
        
        # For 90FV guess
        if guess == "90FV":
            # both letters are incorrect
            if 'F' == pwd[3] or 'V' == pwd[3]:
                return False
                
        # For 68HM guess
        if guess == "68HM":
            # both letters are incorrect
            if 'H' == pwd[3] or 'M' == pwd[3]:
                return False
            # The last letter should be after M
            if pwd[3] <= 'M':
                return False
                
        # For 58RV guess
        if guess == "58RV":
            # letters are incorrect and too late
            if pwd[3] >= 'R':
                return False
                
        # For 47QO guess
        if guess == "47QO":
            # both letters incorrect
            if 'Q' == pwd[3] or 'O' == pwd[3]:
                return False
                
        # For 40SV guess
        if guess == "40SV":
            # letters too late in alphabet
            if pwd[3] >= 'S':
                return False
                
        # For specific letter sequence checking from other guesses
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
            
        # For 57UY guess
        if guess == "57UY":
            # letters are incorrect and too late
            if pwd[3] >= 'P':  # Critical constraint: Y is too late, so our letter must be before P
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