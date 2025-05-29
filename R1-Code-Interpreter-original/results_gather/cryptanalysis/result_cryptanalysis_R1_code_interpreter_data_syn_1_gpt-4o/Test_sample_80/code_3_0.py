# Define the conditions
conditions = [
    ("37ZI", (1, 1, 0, 0)),
    ("74VH", (1, 0, 0, 0)),
    ("39SX", (0, 0, 1, 0)),
    ("67GU", (1, 1, 0, 0)),
    ("86AP", (0, 0, 1, 0)),
    ("80LC", (0, 0, 0, 0)),
    ("38XZ", (0, 0, 0, 1)),
    ("96WO", (0, 0, 0, 0)),
    ("87BK", (1, 1, 0, 0)),
    ("84WL", (0, 0, 0, 0)),
    ("72NV", (2, 0, 0, 0))
]

# Possible numbers and letters
numbers = ['7', '2']
letters = ['S', 'A']

# Check the combination
combination = "72SA"

# Function to check the feedback
def check_feedback(guess, feedback):
    correct_numbers = sum(1 for i in range(2) if guess[i] == combination[i])
    correct_letters = sum(1 for i in range(2, 4) if guess[i] == combination[i])
    wrong_position_numbers = sum(1 for i in range(2) if guess[i] in combination[:2] and guess[i] != combination[i])
    wrong_position_letters = sum(1 for i in range(2, 4) if guess[i] in combination[2:] and guess[i] != combination[i])
    return (correct_numbers + wrong_position_numbers, correct_letters + wrong_position_letters, correct_numbers, correct_letters)

# Verify all conditions
all_conditions_satisfied = all(check_feedback(guess, feedback) == feedback for guess, feedback in conditions)

print(all_conditions_satisfied)