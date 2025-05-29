from itertools import permutations

# Possible numbers and letters
possible_numbers = ['0', '4', '6', '7']
possible_letters = ['N', 'U']

# Generate all permutations of two numbers and two letters
number_permutations = permutations(possible_numbers, 2)
letter_permutations = permutations(possible_letters, 2)

# Check each combination against the conditions
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        # Construct the guess
        guess = num_perm + let_perm
        
        # Check conditions
        if (guess[0] in ['6', '7'] and guess[1] in ['0', '4'] and
            guess[2] == 'N' and guess[3] == 'U'):
            print(guess)
            break