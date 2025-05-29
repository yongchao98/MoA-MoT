from itertools import permutations

# Possible numbers and letters based on deductions
possible_numbers = ['6', '7']
possible_letters = ['G', 'H']

# Generate all permutations of the numbers and letters
number_permutations = list(permutations(possible_numbers, 2))
letter_permutations = list(permutations(possible_letters, 2))

# Check each combination against the feedback
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        combination = num_perm + let_perm
        # Check against the feedback
        if (('7' in combination and '6' in combination) and
            ('G' in combination and 'H' in combination)):
            print(combination)