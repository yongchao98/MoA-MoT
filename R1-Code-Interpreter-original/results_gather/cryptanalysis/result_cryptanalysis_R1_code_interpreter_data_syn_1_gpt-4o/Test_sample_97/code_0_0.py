# Define the possible numbers and letters based on deductions
possible_numbers = [9, 7]
possible_letters = ['T', 'W']

# Function to check the conditions
def check_conditions(numbers, letters):
    conditions = [
        (71, 'BQ', (1, 1, 0, 0)),
        (29, 'BQ', (1, 1, 0, 0)),
        (69, 'BE', (1, 1, 0, 0)),
        (40, 'QT', (0, 0, 1, 1)),
        (48, 'NP', (0, 0, 0, 0)),
        (15, 'FQ', (0, 0, 0, 0)),
        (83, 'PF', (0, 0, 0, 0)),
        (96, 'MJ', (1, 1, 0, 0)),
        (32, 'YJ', (0, 0, 0, 0)),
        (16, 'SU', (0, 0, 0, 0)),
        (89, 'XJ', (1, 1, 0, 0)),
        (36, 'XJ', (0, 0, 0, 0)),
        (35, 'ZU', (0, 0, 0, 0)),
        (36, 'FV', (0, 0, 0, 0))
    ]
    
    for condition in conditions:
        guess_numbers, guess_letters, feedback = condition
        guess_numbers = [int(d) for d in str(guess_numbers)]
        guess_letters = list(guess_letters)
        
        # Check numbers
        correct_numbers = sum(1 for n in guess_numbers if n in numbers)
        correct_positions_numbers = sum(1 for i, n in enumerate(guess_numbers) if n == numbers[i])
        
        # Check letters
        correct_letters = sum(1 for l in guess_letters if l in letters)
        correct_positions_letters = sum(1 for i, l in enumerate(guess_letters) if l == letters[i])
        
        if (correct_numbers, correct_positions_numbers, correct_letters, correct_positions_letters) != feedback:
            return False
    return True

# Check the conditions with the deduced numbers and letters
result = check_conditions(possible_numbers, possible_letters)
print(result)