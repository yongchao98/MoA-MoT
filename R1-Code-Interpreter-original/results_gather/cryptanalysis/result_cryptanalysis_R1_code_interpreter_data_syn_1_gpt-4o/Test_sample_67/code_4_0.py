# Define the correct combination
correct_combination = (4, 5, 'L', 'M')

# Define the feedback for each guess
guesses = [
    (7, 9, 'F', 'V', False),
    (3, 2, 'P', 'Z', False),
    (0, 9, 'E', 'F', False),
    (5, 8, 'Q', 'D', False),
    (7, 9, 'O', 'J', False),
    (6, 4, 'L', 'Y', False),
    (4, 8, 'H', 'I', False),
    (4, 5, 'T', 'G', True),
    (3, 1, 'I', 'B', False),
    (9, 4, 'V', 'W', False),
    (7, 0, 'X', 'N', False),
    (7, 0, 'B', 'I', False),
    (8, 9, 'U', 'G', False),
    (7, 0, 'K', 'G', False)
]

# Function to check if a guess matches the feedback
def check_guess(guess, correct_combination, should_match):
    numbers_match = (guess[0] == correct_combination[0] and guess[1] == correct_combination[1])
    letters_match = (guess[2] == correct_combination[2] and guess[3] == correct_combination[3])
    return (numbers_match and letters_match) == should_match

# Check all guesses
all_conditions_met = all(check_guess(guess, correct_combination, should_match) for *guess, should_match in guesses)

print(all_conditions_met)