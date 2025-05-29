from itertools import permutations

# Possible numbers and letters
numbers = ['1', '6']
letters = ['D', 'F']

# Generate all permutations of the numbers and letters
number_permutations = list(permutations(numbers, 2))
letter_permutations = list(permutations(letters, 2))

# Check each combination
for num_perm in number_permutations:
    for let_perm in letter_permutations:
        combination = num_perm + let_perm
        # Check against all conditions
        if (combination[0] == '1' and combination[1] != '8' and
            combination[2] == 'F' and combination[3] == 'D'):
            print(combination)