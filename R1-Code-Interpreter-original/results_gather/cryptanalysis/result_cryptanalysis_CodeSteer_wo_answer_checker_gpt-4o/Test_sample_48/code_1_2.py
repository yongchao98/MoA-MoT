# Initialize possible values for each position
possible_numbers = [set(range(10)), set(range(10))]
possible_letters = [set('ABCDEFGHIJKLMNOPQRSTUVWXYZ'), set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')]

# Process each guess and feedback
guesses = [
    ("57RU", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("26ZE", "one number is correct but in the wrong position; one number is incorrect and too small; one letter is correct and in the correct position; one letter is incorrect and too early in the alphabet."),
    ("01TK", "both numbers are incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("56GR", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet."),
    ("25XY", "both numbers are incorrect and too small; one letter is correct but in the wrong position; one letter is incorrect and incorrect."),
    ("94KA", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect and too early in the alphabet.")
]

for guess, feedback in guesses:
    numbers = [int(guess[0]), int(guess[1])]
    letters = [guess[2], guess[3]]

    # Process number feedback
    if "both numbers are incorrect" in feedback:
        for i in range(2):
            possible_numbers[i].discard(numbers[i])
    if "one number is correct but in the wrong position" in feedback:
        for i in range(2):
            if numbers[i] in possible_numbers[i]:
                possible_numbers[1-i].add(numbers[i])
                possible_numbers[i].discard(numbers[i])
    if "one number is incorrect and too small" in feedback:
        for i in range(2):
            if numbers[i] in possible_numbers[i]:
                possible_numbers[i] = {n for n in possible_numbers[i] if n > numbers[i]}
    if "both numbers are incorrect and too small" in feedback:
        for i in range(2):
            possible_numbers[i] = {n for n in possible_numbers[i] if n > numbers[i]}

    # Process letter feedback
    if "both letters are incorrect and too early in the alphabet" in feedback:
        for i in range(2):
            possible_letters[i] = {l for l in possible_letters[i] if l > letters[i]}
    if "one letter is correct and in the correct position" in feedback:
        for i in range(2):
            if letters[i] in possible_letters[i]:
                possible_letters[i] = {letters[i]}
    if "one letter is correct but in the wrong position" in feedback:
        for i in range(2):
            if letters[i] in possible_letters[i]:
                possible_letters[1-i].add(letters[i])
                possible_letters[i].discard(letters[i])

# Deduce the final password
final_numbers = [next(iter(possible_numbers[i])) for i in range(2)]
final_letters = [next(iter(possible_letters[i])) for i in range(2)]

# Output the deduced password
password = [str(final_numbers[0]), str(final_numbers[1]), final_letters[0], final_letters[1]]
print(password)