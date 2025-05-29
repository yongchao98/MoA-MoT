from itertools import permutations

# Possible numbers and letters
numbers = ['0', '8']
letters = ['H', 'Q']

# Generate all permutations of the numbers and letters
number_permutations = list(permutations(numbers, 2))
letter_permutations = list(permutations(letters, 2))

# Check each combination
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        combination = list(num_perm) + list(let_perm)
        # Check conditions
        if (combination[0] == '0' and combination[1] == '8' and combination[2] == 'H' and combination[3] == 'Q'):
            print(combination)