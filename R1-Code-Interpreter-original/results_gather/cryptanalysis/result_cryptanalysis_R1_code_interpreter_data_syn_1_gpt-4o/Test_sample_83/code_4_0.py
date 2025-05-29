# Define the conditions
conditions = [
    {"guess": "02PG", "feedback": "both numbers are incorrect; both letters are incorrect and too early in the alphabet."},
    {"guess": "50RS", "feedback": "both numbers are incorrect; both letters are incorrect and too early in the alphabet."},
    {"guess": "76SU", "feedback": "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet."},
    {"guess": "73ZJ", "feedback": "both numbers are incorrect; both letters are incorrect."},
    {"guess": "69WF", "feedback": "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect."},
    {"guess": "37AQ", "feedback": "both numbers are incorrect; both letters are incorrect and too early in the alphabet."},
    {"guess": "05DO", "feedback": "both numbers are incorrect; both letters are incorrect and too early in the alphabet."},
    {"guess": "68PI", "feedback": "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet."},
    {"guess": "98VO", "feedback": "both numbers are incorrect and too large; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."},
    {"guess": "60YG", "feedback": "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect."}
]

# Define the correct password
correct_password = "64OX"

# Function to check if a guess satisfies the feedback
def check_guess(guess, feedback, correct_password):
    numbers_correct = sum(1 for i in range(2) if guess[i] == correct_password[i])
    letters_correct = sum(1 for i in range(2, 4) if guess[i] == correct_password[i])
    
    if "both numbers are incorrect" in feedback:
        if numbers_correct != 0:
            return False
    if "one number is correct and in the correct position" in feedback:
        if numbers_correct != 1:
            return False
    if "one number is correct but in the wrong position" in feedback:
        if numbers_correct != 0:
            return False
    if "both letters are incorrect" in feedback:
        if letters_correct != 0:
            return False
    if "one letter is correct and in the correct position" in feedback:
        if letters_correct != 1:
            return False
    return True

# Check all conditions
all_conditions_satisfied = all(check_guess(cond["guess"], cond["feedback"], correct_password) for cond in conditions)

print(all_conditions_satisfied)