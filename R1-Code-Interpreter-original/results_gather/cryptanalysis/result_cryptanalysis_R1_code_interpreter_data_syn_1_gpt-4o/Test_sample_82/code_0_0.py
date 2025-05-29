from itertools import permutations

# Possible numbers and letters based on deductions
possible_numbers = ['0', '1', '6']
possible_letters = ['M', 'P']

# Generate all permutations of the numbers and letters
number_permutations = list(permutations(possible_numbers, 2))
letter_permutations = list(permutations(possible_letters, 2))

# Check each combination against the conditions
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        # Construct the guess
        guess = num_perm + let_perm
        
        # Check against all conditions
        if (guess[0] != '8' and guess[1] != '7' and
            guess[0] != '7' and guess[1] != '1' and
            guess[0] != '2' and guess[1] != '3' and
            guess[0] != '5' and guess[1] != '3' and
            guess[0] != '4' and guess[1] != '3' and
            guess[0] != '0' and guess[1] != '1' and
            guess[2] != 'C' and guess[3] != 'F' and
            guess[2] != 'S' and guess[3] != 'M' and
            guess[2] != 'A' and guess[3] != 'Y' and
            guess[2] != 'P' and guess[3] != 'D' and
            guess[2] != 'U' and guess[3] != 'F' and
            guess[2] != 'U' and guess[3] != 'S'):
            print(guess)