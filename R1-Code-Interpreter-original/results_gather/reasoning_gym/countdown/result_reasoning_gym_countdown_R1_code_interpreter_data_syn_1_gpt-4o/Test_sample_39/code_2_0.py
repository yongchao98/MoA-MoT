# Check if any combination of the remaining numbers can form 50
from itertools import permutations

numbers = [51, 4, 32, 13]
target = 50

# Try all permutations of numbers
for num_perm in permutations(numbers):
    # Check if any combination of operations results in the target
    if num_perm[0] - num_perm[1] + num_perm[2] - num_perm[3] == target:
        print(f"{num_perm[0]} - {num_perm[1]} + {num_perm[2]} - {num_perm[3]}")
        break