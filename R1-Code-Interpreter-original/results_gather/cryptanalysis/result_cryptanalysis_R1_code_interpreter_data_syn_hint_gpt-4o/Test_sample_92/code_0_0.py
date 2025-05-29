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

# Deduce the correct password
correct_numbers = [2, 6]
correct_letters = ['X', 'V']

# Verify the deductions
for guess, feedback in guesses:
    numbers = [int(guess[0]), int(guess[1])]
    letters = [guess[2], guess[3]]
    
    # Check numbers
    if "both numbers are incorrect" in feedback:
        assert all(num not in correct_numbers for num in numbers)
    if "one number is correct and in the correct position" in feedback:
        assert any(num == correct_numbers[i] for i, num in enumerate(numbers))
    if "one number is correct but in the wrong position" in feedback:
        assert any(num in correct_numbers and num != correct_numbers[i] for i, num in enumerate(numbers))
    
    # Check letters
    if "both letters are incorrect" in feedback:
        assert all(letter not in correct_letters for letter in letters)
    if "one letter is correct but in the wrong position" in feedback:
        assert any(letter in correct_letters and letter != correct_letters[i] for i, letter in enumerate(letters))

print("The correct password is:", correct_numbers + correct_letters)