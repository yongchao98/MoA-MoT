from itertools import permutations

# Possible numbers and letters
numbers = [0, 1, 2, 3]
letters = ['I', 'M']

# Generate all permutations of two numbers and two letters
possible_combinations = permutations(numbers, 2)
possible_letters = permutations(letters, 2)

# Check each combination against the conditions
for num_comb in possible_combinations:
    for let_comb in possible_letters:
        # Construct the guess
        guess = list(num_comb) + list(let_comb)
        
        # Check conditions
        if (guess[0] != 0 and guess[1] == 0) and (guess[2] == 'M' and guess[3] == 'I'):
            print(guess)
            break