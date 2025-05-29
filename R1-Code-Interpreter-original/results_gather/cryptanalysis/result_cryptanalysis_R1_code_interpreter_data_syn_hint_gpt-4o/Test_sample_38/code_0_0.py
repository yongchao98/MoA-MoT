from itertools import permutations

# Possible numbers and letters based on the conditions
possible_numbers = [0, 2, 3]
possible_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

# Check each permutation of numbers and letters
for num_perm in permutations(possible_numbers, 2):
    for letter_perm in permutations(possible_letters, 2):
        # Construct the guess
        guess = list(num_perm) + list(letter_perm)
        
        # Check against the conditions
        if (guess[0] != 4 and guess[1] != 7 and
            guess[0] != 9 and guess[1] != 1 and
            (guess[2] == 'L' or guess[3] == 'L') and
            (guess[2] == 'K' or guess[3] == 'K')):
            print(guess)
            break