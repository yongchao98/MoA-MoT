# Define the guesses and feedback
guesses = [
    ("07UC", (1, 1, 0)),  # one number correct but wrong position, one number incorrect and too large, both letters incorrect
    ("60FG", (1, 1, 0)),  # one number correct and correct position, one number incorrect and too large, both letters incorrect
    ("59VF", (0, 0, 1)),  # both numbers incorrect and too large, one letter correct and correct position, one letter incorrect
    ("85EK", (0, 0, 1)),  # both numbers incorrect and too large, one letter correct but wrong position, one letter incorrect
    ("18MA", (1, 1, 0))   # one number correct and correct position, one number incorrect and too large, both letters incorrect
]

# Possible numbers and letters
numbers = "0123456789"
letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Function to check if a guess matches the feedback
def check_guess(guess, feedback, password):
    num_correct_pos = sum(guess[i] == password[i] for i in range(2))
    num_correct_wrong_pos = sum(guess[i] in password[:2] and guess[i] != password[i] for i in range(2))
    letter_correct_pos = sum(guess[i] == password[i] for i in range(2, 4))
    letter_correct_wrong_pos = sum(guess[i] in password[2:] and guess[i] != password[i] for i in range(2, 4))
    
    return (num_correct_pos + num_correct_wrong_pos == feedback[0] and
            num_correct_pos == feedback[1] and
            letter_correct_pos + letter_correct_wrong_pos == feedback[2])

# Iterate over all possible combinations
for n1 in numbers:
    for n2 in numbers:
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                password = n1 + n2 + l1 + l2
                if all(check_guess(guess, feedback, password) for guess, feedback in guesses):
                    print(password)