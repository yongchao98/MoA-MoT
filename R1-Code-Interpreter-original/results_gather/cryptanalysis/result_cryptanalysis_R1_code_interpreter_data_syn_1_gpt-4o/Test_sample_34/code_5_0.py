# Define the possible password based on deductions
password = ["3", "6", "V", "C"]

# Define the guesses and feedback
guesses = [
    ("18VC", "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."),
    ("54OF", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("36AV", "both numbers are correct and in the correct positions; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet."),
    ("42GE", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("08NQ", "both numbers are incorrect; both letters are incorrect.")
]

# Function to check if the password satisfies all conditions
def check_password(password, guesses):
    for guess, feedback in guesses:
        numbers_correct = (password[0] == guess[0] and password[1] == guess[1])
        letters_correct = (password[2] == guess[2] and password[3] == guess[3])
        
        if "both numbers are incorrect" in feedback and numbers_correct:
            return False
        if "both numbers are correct" in feedback and not numbers_correct:
            return False
        if "one letter is correct and in the correct position" in feedback:
            if not ((password[2] == guess[2] and password[3] != guess[3]) or (password[3] == guess[3] and password[2] != guess[2])):
                return False
        if "one letter is correct but in the wrong position" in feedback:
            if not ((password[2] == guess[3] or password[3] == guess[2]) and not letters_correct):
                return False
        if "both letters are incorrect" in feedback:
            if password[2] in guess or password[3] in guess:
                return False
        if "both letters are incorrect and too early in the alphabet" in feedback:
            if password[2] <= guess[2] or password[3] <= guess[3]:
                return False
    return True

# Check if the password satisfies all conditions
is_valid = check_password(password, guesses)
print(is_valid)