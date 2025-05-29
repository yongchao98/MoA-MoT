# Initialize possible numbers and letters
possible_numbers = set(range(10))
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Initialize known positions
number_positions = [None, None]
letter_positions = [None, None]

# Define the guesses and feedback
guesses = [
    ("52OL", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("20FP", "one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet"),
    ("08EX", "one number is correct but in the wrong position; one number is incorrect and too large; both letters are incorrect"),
    ("58RN", "both numbers are incorrect and too large; both letters are incorrect and too early in the alphabet"),
    ("64UM", "one number is correct but in the wrong position; one number is incorrect and too large; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet"),
    ("02PV", "one number is correct but in the wrong position; one number is incorrect and incorrect; both letters are incorrect"),
    ("62CO", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("91CD", "both numbers are incorrect; both letters are incorrect and too early in the alphabet"),
    ("96SR", "both numbers are incorrect and too large; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet")
]

# Process each guess
for guess, feedback in guesses:
    numbers = [int(guess[0]), int(guess[1])]
    letters = [guess[2], guess[3]]

    if "both numbers are incorrect" in feedback:
        possible_numbers -= set(numbers)
    if "both letters are incorrect" in feedback:
        possible_letters -= set(letters)
    if "one number is correct and in the correct position" in feedback:
        for i, num in enumerate(numbers):
            if num in possible_numbers:
                number_positions[i] = num
    if "one number is correct but in the wrong position" in feedback:
        for i, num in enumerate(numbers):
            if num in possible_numbers:
                number_positions[1 - i] = num
    if "one letter is correct and in the correct position" in feedback:
        for i, letter in enumerate(letters):
            if letter in possible_letters:
                letter_positions[i] = letter
    if "one letter is correct but in the wrong position" in feedback:
        for i, letter in enumerate(letters):
            if letter in possible_letters:
                letter_positions[1 - i] = letter

# Deduce the final password
password = [
    str(number_positions[0] if number_positions[0] is not None else next(iter(possible_numbers))),
    str(number_positions[1] if number_positions[1] is not None else next(iter(possible_numbers))),
    letter_positions[0] if letter_positions[0] is not None else next(iter(possible_letters)),
    letter_positions[1] if letter_positions[1] is not None else next(iter(possible_letters))
]

print(f"<<< {password} >>>")