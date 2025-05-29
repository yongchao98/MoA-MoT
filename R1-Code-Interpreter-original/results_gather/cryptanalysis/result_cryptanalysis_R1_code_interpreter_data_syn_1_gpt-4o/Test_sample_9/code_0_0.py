from itertools import permutations

# Possible numbers and letters
possible_numbers = [0, 3, 7, 8]
possible_letters = ['E', 'W']

# Find the correct combination
for num_perm in permutations(possible_numbers, 2):
    if 7 in num_perm:
        for letter_perm in permutations(possible_letters, 2):
            # Check if the number 7 is in the wrong position
            if num_perm[0] != 7:
                # Print the correct combination
                print([num_perm[0], num_perm[1], letter_perm[0], letter_perm[1]])