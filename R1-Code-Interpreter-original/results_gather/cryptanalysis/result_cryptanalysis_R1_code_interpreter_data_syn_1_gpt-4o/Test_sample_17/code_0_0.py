from itertools import permutations

# Possible numbers and letters based on the analysis
numbers = [1, 6]
letters = ['Q', 'X', 'R', 'S', 'T', 'U', 'V', 'W', 'Z']

# Generate all permutations of the numbers and letters
number_permutations = permutations(numbers)
letter_permutations = permutations(letters, 2)

# Check each combination
for num_perm in number_permutations:
    for letter_perm in letter_permutations:
        # Construct the guess
        guess = list(num_perm) + list(letter_perm)
        
        # Check the conditions
        if (guess[0] == 6 and guess[1] == 1 and
            ('Q' in guess[2:] or 'X' in guess[2:])):
            print(guess)
            break