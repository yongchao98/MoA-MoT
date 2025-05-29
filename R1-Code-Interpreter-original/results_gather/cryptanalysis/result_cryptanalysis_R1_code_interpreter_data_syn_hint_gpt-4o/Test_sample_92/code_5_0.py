# Define the guesses and feedback
guesses = [
    ("58TD", "both numbers are incorrect; both letters are incorrect."),
    ("52SK", "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect."),
    ("19UI", "both numbers are incorrect; both letters are incorrect."),
    ("75AC", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("86WD", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect."),
    ("72KL", "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet."),
    ("08LP", "both numbers are incorrect; both letters are incorrect."),
    ("51PB", "both numbers are incorrect; both letters are incorrect."),
    ("96FX", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect."),
    ("02HU", "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect."),
    ("96FJ", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect and too early in the alphabet."),
    ("12KC", "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("92TG", "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect."),
    ("98OJ", "both numbers are incorrect and too large; both letters are incorrect."),
    ("06ML", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("02VX", "one number is correct and in the correct position; one number is incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and too late in the alphabet.")
]

# Possible numbers and letters
possible_numbers = [str(i) for i in range(10)]
possible_letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]

# Function to check if a guess matches the feedback
def matches_feedback(guess, feedback, correct_numbers, correct_letters):
    numbers = [guess[0], guess[1]]
    letters = [guess[2], guess[3]]
    
    # Check numbers
    if "both numbers are incorrect" in feedback:
        if any(num in correct_numbers for num in numbers):
            return False
    if "one number is correct and in the correct position" in feedback:
        if not any(num == correct_numbers[i] for i, num in enumerate(numbers)):
            return False
    if "one number is correct but in the wrong position" in feedback:
        if not any(num in correct_numbers and num != correct_numbers[i] for i, num in enumerate(numbers)):
            return False
    
    # Check letters
    if "both letters are incorrect" in feedback:
        if any(letter in correct_letters for letter in letters):
            return False
    if "one letter is correct but in the wrong position" in feedback:
        if not any(letter in correct_letters and letter != correct_letters[i] for i, letter in enumerate(letters)):
            return False
    
    return True

# Iterate through all possible combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 == num2:
            continue
        for letter1 in possible_letters:
            for letter2 in possible_letters:
                if letter1 == letter2:
                    continue
                correct_numbers = [num1, num2]
                correct_letters = [letter1, letter2]
                if all(matches_feedback(guess, feedback, correct_numbers, correct_letters) for guess, feedback in guesses):
                    print("The correct password is:", correct_numbers + correct_letters)
                    break