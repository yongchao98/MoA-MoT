# Define the guesses and feedback
guesses = [
    ("92BJ", (0, 0)),
    ("06RU", (1, 0)),
    ("84CV", (0, 1)),
    ("75VY", (1, 0)),
    ("16UD", (1, 0)),
    ("98NG", (0, 0)),
    ("03QD", (0, 1))
]

# Define the possible numbers and letters
possible_numbers = set("0123456789")
possible_letters = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

# Eliminate impossible numbers and letters based on feedback
for guess, feedback in guesses:
    numbers, letters = guess[:2], guess[2:]
    num_correct, let_correct = feedback
    
    if num_correct == 0:
        possible_numbers -= set(numbers)
    if let_correct == 0:
        possible_letters -= set(letters)

# Further refine based on specific feedback
# From guess 6, numbers must be less than 8
possible_numbers -= set("89")

# From guess 2, one number is correct and in the correct position
# Since 0 is too small, 6 must be correct
# From guess 4, one number is correct and in the correct position
# Since 5 is too small, 7 must be correct
# From guess 5, one number is correct and in the correct position
# Since 1 is too small, 6 must be correct
correct_numbers = "67"

# From guess 3, one letter is correct and in the correct position
# Since V is incorrect, C must be correct
# From guess 7, one letter is correct but in the wrong position
# Since D is incorrect, Q must be correct
correct_letters = "CQ"

# Verify the solution
solution = "67CQ"
is_valid = True
for guess, feedback in guesses:
    num_correct, let_correct = feedback
    numbers, letters = guess[:2], guess[2:]
    
    # Check numbers
    num_match = sum(1 for i in range(2) if solution[i] == numbers[i])
    # Check letters
    let_match = sum(1 for i in range(2) if solution[i+2] == letters[i])
    
    if num_match != num_correct or let_match != let_correct:
        is_valid = False
        break

print("Solution is valid:", is_valid)
print("Password:", list(solution))