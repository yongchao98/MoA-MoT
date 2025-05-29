# Define the guesses and feedback
guesses = [
    ("20BH", (1, 1, 0)),  # (correct numbers in wrong position, incorrect numbers, incorrect letters)
    ("85IM", (0, 2, 0)),  # (correct numbers, incorrect numbers, incorrect letters)
    ("09OQ", (2, 0, 1)),  # (correct numbers in correct position, incorrect numbers, correct letters in wrong position)
    ("28BR", (0, 2, 1))   # (correct numbers, incorrect numbers, correct letters in wrong position)
]

# Possible numbers and letters
numbers = "0123456789"
letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Deduce the correct combination
correct_numbers = "09"
correct_letters = "RO"

# Check the combination against the feedback
def check_combination(guess, feedback):
    num_correct_pos = sum(1 for i in range(2) if guess[i] == correct_numbers[i])
    num_correct_wrong_pos = sum(1 for i in range(2) if guess[i] in correct_numbers and guess[i] != correct_numbers[i])
    num_incorrect = 2 - num_correct_pos - num_correct_wrong_pos
    
    letter_correct_pos = sum(1 for i in range(2, 4) if guess[i] == correct_letters[i-2])
    letter_correct_wrong_pos = sum(1 for i in range(2, 4) if guess[i] in correct_letters and guess[i] != correct_letters[i-2])
    letter_incorrect = 2 - letter_correct_pos - letter_correct_wrong_pos
    
    return (num_correct_wrong_pos, num_incorrect, letter_correct_wrong_pos) == feedback

# Validate the solution
valid = all(check_combination(guess, feedback) for guess, feedback in guesses)

print(valid)