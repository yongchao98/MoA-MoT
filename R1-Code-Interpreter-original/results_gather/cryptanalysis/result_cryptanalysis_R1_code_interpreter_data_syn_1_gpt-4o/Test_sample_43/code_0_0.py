# Define the guesses and feedback
guesses = [
    ("61OJ", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("98QN", "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet."),
    ("58FC", "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet."),
    ("72YK", "both numbers are incorrect; both letters are incorrect."),
    ("02ER", "both numbers are incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("35QW", "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect."),
    ("71KJ", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("79PX", "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."),
    ("86TV", "one number is correct but in the wrong position; one number is incorrect and incorrect; one letter is correct and in the correct position; one letter is incorrect and incorrect.")
]

# Define the correct combination
correct_combination = "83TX"

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback):
    numbers_correct = 0
    numbers_wrong_position = 0
    letters_correct = 0
    letters_wrong_position = 0
    
    # Check numbers
    for i in range(2):
        if guess[i] == correct_combination[i]:
            numbers_correct += 1
        elif guess[i] in correct_combination[:2]:
            numbers_wrong_position += 1
    
    # Check letters
    for i in range(2, 4):
        if guess[i] == correct_combination[i]:
            letters_correct += 1
        elif guess[i] in correct_combination[2:]:
            letters_wrong_position += 1
    
    # Verify feedback
    if "both numbers are incorrect" in feedback:
        if numbers_correct != 0:
            return False
    if "one number is correct and in the correct position" in feedback:
        if numbers_correct != 1:
            return False
    if "one number is correct but in the wrong position" in feedback:
        if numbers_wrong_position != 1:
            return False
    if "both letters are incorrect" in feedback:
        if letters_correct != 0:
            return False
    if "one letter is correct and in the correct position" in feedback:
        if letters_correct != 1:
            return False
    if "one letter is correct but in the wrong position" in feedback:
        if letters_wrong_position != 1:
            return False
    
    return True

# Check all guesses
all_correct = all(check_guess(guess, feedback) for guess, feedback in guesses)

print(all_correct)