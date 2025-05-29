# Define the guesses and feedback
guesses = [
    ("74SE", "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet"),
    ("76OX", "both numbers are incorrect; both letters are incorrect"),
    ("94LJ", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("28HI", "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet"),
    ("27LJ", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("36WF", "both numbers are incorrect; both letters are incorrect"),
    ("71IS", "both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet"),
    ("12YB", "both numbers are incorrect and too small; both letters are incorrect"),
    ("87WE", "one number is correct but in the wrong position; one number is incorrect and incorrect; both letters are incorrect"),
    ("25YP", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect"),
    ("84UL", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect"),
    ("30VF", "both numbers are incorrect and too small; both letters are incorrect"),
    ("68QG", "one number is correct and in the correct position; one number is incorrect and incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet")
]

# Define the correct password
correct_password = ["8", "9", "S", "E"]

# Function to check if a guess matches the feedback
def check_guess(guess, feedback):
    numbers = guess[:2]
    letters = guess[2:]
    
    # Check numbers
    if "both numbers are incorrect" in feedback:
        if any(n in correct_password[:2] for n in numbers):
            return False
    if "one number is correct and in the correct position" in feedback:
        if not any(correct_password[i] == numbers[i] for i in range(2)):
            return False
    if "one number is correct but in the wrong position" in feedback:
        if not any(correct_password[i] != numbers[i] and numbers[i] in correct_password[:2] for i in range(2)):
            return False
    if "one number is incorrect and too small" in feedback:
        if not any(int(n) < int(correct_password[i]) for i, n in enumerate(numbers)):
            return False
    
    # Check letters
    if "both letters are incorrect" in feedback:
        if any(l in correct_password[2:] for l in letters):
            return False
    if "one letter is correct and in the correct position" in feedback:
        if not any(correct_password[i+2] == letters[i] for i in range(2)):
            return False
    if "one letter is correct but in the wrong position" in feedback:
        if not any(correct_password[i+2] != letters[i] and letters[i] in correct_password[2:] for i in range(2)):
            return False
    if "one letter is incorrect and too early in the alphabet" in feedback:
        if not any(ord(l) < ord(correct_password[i+2]) for i, l in enumerate(letters)):
            return False
    
    return True

# Check all guesses
all_correct = all(check_guess(guess, feedback) for guess, feedback in guesses)

# Print the result
print(all_correct)