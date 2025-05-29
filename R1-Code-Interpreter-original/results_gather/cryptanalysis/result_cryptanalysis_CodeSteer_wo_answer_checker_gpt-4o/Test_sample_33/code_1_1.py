from itertools import permutations

def simulate_feedback(guess, combination):
    feedback = {'correct_position': 0, 'incorrect_position': 0, 'too_small': 0, 'incorrect': 0}
    
    # Check numbers
    for i in range(2):
        if guess[i] == combination[i]:
            feedback['correct_position'] += 1
        elif guess[i] in combination[:2]:
            feedback['incorrect_position'] += 1
        elif guess[i] < combination[i]:
            feedback['too_small'] += 1
        else:
            feedback['incorrect'] += 1
    
    # Check letters
    for i in range(2, 4):
        if guess[i] == combination[i]:
            feedback['correct_position'] += 1
        elif guess[i] in combination[2:]:
            feedback['incorrect_position'] += 1
        else:
            feedback['incorrect'] += 1
    
    return feedback

def is_valid_combination(combination, guesses):
    for guess, expected_feedback in guesses:
        feedback = simulate_feedback(guess, combination)
        if feedback != expected_feedback:
            return False
    return True

# Define the guesses and expected feedback
guesses = [
    ("35WJ", {'correct_position': 1, 'too_small': 1, 'incorrect': 2}),
    ("67VX", {'incorrect': 4}),
    ("37QK", {'incorrect': 2, 'correct_position': 1, 'incorrect_position': 1}),
    ("75GT", {'correct_position': 1, 'incorrect': 3}),
    ("34QE", {'too_small': 2, 'incorrect': 2}),
    ("27CU", {'incorrect': 4}),
    ("56FA", {'incorrect_position': 1, 'incorrect': 3}),
    ("84YT", {'correct_position': 2, 'too_small': 1, 'incorrect': 1})
]

# Generate all possible combinations of two numbers and two letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
possible_combinations = [''.join(num_perm) + ''.join(let_perm) 
                         for num_perm in permutations(numbers, 2) 
                         for let_perm in permutations(letters, 2)]

# Find the valid combination
for combination in possible_combinations:
    if is_valid_combination(combination, guesses):
        print(f"<<< {list(combination)} >>>")
        break