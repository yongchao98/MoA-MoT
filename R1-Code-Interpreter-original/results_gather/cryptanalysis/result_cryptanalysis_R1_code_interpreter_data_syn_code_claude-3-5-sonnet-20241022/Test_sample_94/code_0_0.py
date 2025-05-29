from itertools import product
import string

def is_valid_password(password, guesses_feedback):
    # Helper function to check if a number is "too large" or "too small"
    def compare_numbers(guess_num, actual_num, feedback):
        if "too large" in feedback:
            return int(guess_num) > int(actual_num)
        if "too small" in feedback:
            return int(guess_num) < int(actual_num)
        return True

    for guess, feedback in guesses_feedback:
        # Check each position
        correct_count = 0
        wrong_pos_count = 0
        
        # Special handling for guess 61VG which has specific feedback
        if guess == "61VG":
            # One number correct and in position
            if not ((password[0] == '6' and password[1] != '1') or 
                   (password[1] == '1' and password[0] != '6')):
                return False
            # One number incorrect and too small
            if not (int('6' if password[0] != '6' else '1') < 
                   int(password[0] if password[0] != '6' else password[1])):
                return False
            # One letter correct but wrong position
            if not (('V' in password[2:] and password.index('V') != 2) or 
                   ('G' in password[2:] and password.index('G') != 3)):
                return False
            continue

        # Check numbers
        for i in range(2):
            if "both numbers are incorrect" in feedback:
                if guess[i] == password[i]:
                    return False
            elif "one number is correct but in the wrong position" in feedback:
                if guess[i] == password[i]:
                    return False
                if guess[i] in password[:2]:
                    wrong_pos_count += 1
            elif "one number is correct and in the correct position" in feedback:
                if guess[i] == password[i]:
                    correct_count += 1

        # Verify number conditions
        if "both numbers are incorrect and too large" in feedback:
            if not all(int(guess[i]) > int(password[i]) for i in range(2)):
                return False
        elif "both numbers are incorrect and too small" in feedback:
            if not all(int(guess[i]) < int(password[i]) for i in range(2)):
                return False

        # Check letters
        if "both letters are incorrect and too late in the alphabet" in feedback:
            if not all(guess[i] > password[i] for i in range(2, 4)):
                return False
        elif "both letters are incorrect" in feedback:
            if any(guess[i] == password[i] for i in range(2, 4)):
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

# Try all possible combinations
for n1, n2 in product(numbers, repeat=2):
    if n1 != n2:  # No repeating numbers
        for l1, l2 in product(letters, repeat=2):
            if l1 != l2:  # No repeating letters
                password = n1 + n2 + l1 + l2
                if is_valid_password(password, guesses_feedback):
                    valid_passwords.append([n1, n2, l1, l2])

print(valid_passwords)