# Initialize possible candidates for each position
numbers = [{'0', '9'}, {'0', '9'}]  # Two positions for numbers
letters = [{'O', 'R'}, {'O', 'R'}]  # Two positions for letters

# Process each guess and feedback
guesses = [
    ('20BH', 'one number is correct but in the wrong position; one number is incorrect and incorrect; both letters are incorrect and too early in the alphabet.'),
    ('85IM', 'both numbers are incorrect; both letters are incorrect and too early in the alphabet.'),
    ('09OQ', 'both numbers are correct and in the correct positions; one letter is correct but in the wrong position; one letter is incorrect and incorrect.'),
    ('28BR', 'both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and too early in the alphabet.')
]

# Deduce the password
for guess, feedback in guesses:
    num_guess = guess[:2]
    letter_guess = guess[2:]

    if 'both numbers are correct and in the correct positions' in feedback:
        numbers[0] = {num_guess[0]}
        numbers[1] = {num_guess[1]}
    elif 'one number is correct but in the wrong position' in feedback:
        if num_guess[0] in numbers[1]:
            numbers[1] = {num_guess[0]}
        if num_guess[1] in numbers[0]:
            numbers[0] = {num_guess[1]}
    elif 'both numbers are incorrect' in feedback:
        numbers[0].discard(num_guess[0])
        numbers[1].discard(num_guess[1])

    if 'one letter is correct but in the wrong position' in feedback:
        if letter_guess[0] in letters[1]:
            letters[1] = {letter_guess[0]}
        if letter_guess[1] in letters[0]:
            letters[0] = {letter_guess[1]}
    elif 'both letters are incorrect' in feedback:
        letters[0].discard(letter_guess[0])
        letters[1].discard(letter_guess[1])

# Extract the final password
password = [list(numbers[0])[0], list(numbers[1])[0], list(letters[0])[0], list(letters[1])[0]]
print(f"<<< {password} >>>")