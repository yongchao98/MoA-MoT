# Define the correct numbers and letters
correct_numbers = [3, 8]
correct_letters = ['M', 'L']

# Define the guesses and feedback
guesses = [
    (46, 'DL', "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."),
    (31, 'GM', "one number is correct and in the correct position; one number is incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet."),
    (75, 'WL', "both numbers are incorrect; one letter is correct and in the correct position; one letter is incorrect and too late in the alphabet."),
    (74, 'IA', "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    (4, 'OB', "both numbers are incorrect; both letters are incorrect."),
    (38, 'FP', "both numbers are correct and in the correct positions; both letters are incorrect.")
]

# Check each guess against the correct combination
def check_guess(guess, correct_numbers, correct_letters):
    numbers, letters = guess
    num1, num2 = divmod(numbers, 10)
    letter1, letter2 = letters

    # Check numbers
    num_feedback = []
    if num1 == correct_numbers[0]:
        num_feedback.append("first number correct")
    elif num1 in correct_numbers:
        num_feedback.append("first number correct but wrong position")
    else:
        num_feedback.append("first number incorrect")

    if num2 == correct_numbers[1]:
        num_feedback.append("second number correct")
    elif num2 in correct_numbers:
        num_feedback.append("second number correct but wrong position")
    else:
        num_feedback.append("second number incorrect")

    # Check letters
    letter_feedback = []
    if letter1 == correct_letters[0]:
        letter_feedback.append("first letter correct")
    elif letter1 in correct_letters:
        letter_feedback.append("first letter correct but wrong position")
    else:
        letter_feedback.append("first letter incorrect")

    if letter2 == correct_letters[1]:
        letter_feedback.append("second letter correct")
    elif letter2 in correct_letters:
        letter_feedback.append("second letter correct but wrong position")
    else:
        letter_feedback.append("second letter incorrect")

    return num_feedback, letter_feedback

# Validate all guesses
for guess in guesses:
    numbers, letters, feedback = guess
    num_feedback, letter_feedback = check_guess((numbers, letters), correct_numbers, correct_letters)
    print(f"Guess: {numbers}{letters}, Feedback: {feedback}")
    print(f"Computed Feedback: {num_feedback}, {letter_feedback}\n")
