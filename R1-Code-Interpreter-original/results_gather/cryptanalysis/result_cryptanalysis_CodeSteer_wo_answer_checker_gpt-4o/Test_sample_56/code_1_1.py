# Initialize possible values for each position
possible_numbers = [{str(i) for i in range(10)} for _ in range(2)]
possible_letters = [{chr(i) for i in range(ord('A'), ord('Z') + 1)} for _ in range(2)]

# Define the constraints based on the feedback
constraints = [
    ("93YW", "both numbers are correct and in the correct positions; both letters are incorrect and too late in the alphabet."),
    ("34PL", "one number is correct but in the wrong position; one number is incorrect; both letters are incorrect."),
    ("14YH", "both numbers are incorrect; both letters are incorrect."),
    ("70FZ", "both numbers are incorrect; both letters are incorrect."),
    ("46GK", "both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect."),
    ("32RP", "one number is correct but in the wrong position; one number is incorrect and too small; both letters are incorrect."),
    ("63VH", "one number is correct and in the correct position; one number is incorrect; both letters are incorrect."),
    ("13NC", "one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect."),
    ("61BE", "both numbers are incorrect; both letters are incorrect and too early in the alphabet."),
    ("43KU", "one number is correct and in the correct position; one number is incorrect; one letter is correct but in the wrong position; one letter is incorrect.")
]

# Function to apply constraints
def apply_constraints():
    for guess, feedback in constraints:
        numbers, letters = guess[:2], guess[2:]
        
        if "both numbers are correct and in the correct positions" in feedback:
            possible_numbers[0] = {numbers[0]}
            possible_numbers[1] = {numbers[1]}
        
        if "one number is correct but in the wrong position" in feedback:
            if numbers[0] in possible_numbers[1]:
                possible_numbers[1] = {numbers[0]}
            if numbers[1] in possible_numbers[0]:
                possible_numbers[0] = {numbers[1]}
        
        if "both numbers are incorrect" in feedback:
            for i in range(2):
                possible_numbers[i].discard(numbers[0])
                possible_numbers[i].discard(numbers[1])
        
        if "both letters are incorrect" in feedback:
            for i in range(2):
                possible_letters[i].discard(letters[0])
                possible_letters[i].discard(letters[1])
        
        if "one letter is correct but in the wrong position" in feedback:
            if letters[0] in possible_letters[1]:
                possible_letters[1] = {letters[0]}
            if letters[1] in possible_letters[0]:
                possible_letters[0] = {letters[1]}
        
        if "both letters are incorrect and too late in the alphabet" in feedback:
            for i in range(2):
                possible_letters[i] = {l for l in possible_letters[i] if l < letters[0] and l < letters[1]}
        
        if "both letters are incorrect and too early in the alphabet" in feedback:
            for i in range(2):
                possible_letters[i] = {l for l in possible_letters[i] if l > letters[0] and l > letters[1]}

# Apply constraints iteratively
apply_constraints()

# Extract the final password
password = [next(iter(possible_numbers[i])) for i in range(2)] + [next(iter(possible_letters[i])) for i in range(2)]
print(f"<<< {password} >>>")