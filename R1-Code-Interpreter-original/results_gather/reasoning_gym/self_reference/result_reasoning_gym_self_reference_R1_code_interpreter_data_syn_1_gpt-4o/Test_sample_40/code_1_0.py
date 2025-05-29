from itertools import product

# Define the statements
statements = [
    lambda t: t >= 6,  # Statement 1
    lambda f: f <= 5,  # Statement 2
    lambda t: t == 2,  # Statement 3
    lambda f: f == 2,  # Statement 4
    lambda t, f: (t == 2) != (f == 2),  # Statement 5
    lambda t: t in {2, 3, 5, 7},  # Statement 6
    lambda f: f in {4, 6}  # Statement 7
]

# Count the number of valid solutions
valid_solutions = 0

# Iterate over all possible combinations of true (1) and false (0) for 7 statements
for combination in product([0, 1], repeat=7):
    true_count = sum(combination)
    false_count = 7 - true_count
    
    # Check if all statements are satisfied
    if all(statements[i](true_count, false_count) if i in {4} else statements[i](true_count) for i in range(7)):
        valid_solutions += 1

print(valid_solutions)