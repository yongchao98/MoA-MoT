from itertools import permutations

# Possible numbers and letters based on the analysis
possible_numbers = [0, 4, 6, 7]
possible_letters = ['N', 'U']

# Generate all permutations of the numbers and letters
number_permutations = permutations(possible_numbers, 2)
letter_permutations = permutations(possible_letters, 2)

# Check each combination against the conditions
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        # Construct the guess
        guess = list(num_perm) + list(let_perm)
        
        # Check against each condition
        if (guess[0] in [6, 7] and guess[1] in [0, 4] and
            guess[2] == 'N' and guess[3] == 'U'):
            print(guess)
            break