# Initialize possible values
possible_numbers = set(range(10))
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Define the feedback for each guess
guesses = [
    ("57UY", "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    ("40SV", "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    ("01ZD", "both numbers are incorrect and too small; both letters are incorrect."),
    ("58RV", "both numbers are incorrect; both letters are incorrect and too late in the alphabet."),
    ("90FV", "both numbers are incorrect; both letters are incorrect."),
    ("38BY", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect."),
    ("47QO", "both numbers are incorrect; both letters are incorrect."),
    ("94WK", "both numbers are incorrect; both letters are incorrect."),
    ("68HM", "one number is correct and in the correct position; one number is incorrect and too large; both letters are incorrect."),
    ("87HB", "both numbers are incorrect and too large; both letters are incorrect."),
    ("97SH", "both numbers are incorrect and too large; both letters are incorrect."),
    ("29TE", "both numbers are incorrect; both letters are incorrect."),
    ("18HS", "both numbers are incorrect; both letters are incorrect."),
    ("92LU", "both numbers are incorrect; both letters are incorrect."),
    ("87FM", "both numbers are incorrect and too large; both letters are incorrect."),
    ("61VG", "one number is correct and in the correct position; one number is incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and too late in the alphabet.")
]

# Process each guess
for guess, feedback in guesses:
    numbers = set(int(guess[i]) for i in range(2))
    letters = set(guess[i] for i in range(2, 4))
    
    if "both numbers are incorrect" in feedback:
        possible_numbers -= numbers
    if "both letters are incorrect" in feedback:
        possible_letters -= letters
    if "one number is correct but in the wrong position" in feedback:
        possible_numbers &= numbers
    if "one number is correct and in the correct position" in feedback:
        possible_numbers &= numbers
    if "one letter is correct but in the wrong position" in feedback:
        possible_letters &= letters

# Deduce the correct numbers and letters
correct_numbers = list(possible_numbers)
correct_letters = list(possible_letters)

# Output the deduced password
password = [str(correct_numbers[0]), str(correct_numbers[1]), correct_letters[0], correct_letters[1]]
print(f"<<< {password} >>>")