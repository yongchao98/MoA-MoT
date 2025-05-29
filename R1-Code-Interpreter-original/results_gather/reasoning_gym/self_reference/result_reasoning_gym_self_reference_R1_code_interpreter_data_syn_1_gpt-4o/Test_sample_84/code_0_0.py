from itertools import combinations

# Statements
statements = [
    lambda t, f: t >= 3,  # Statement 1
    lambda t, f: f <= 3,  # Statement 2
    lambda t, f: t == 4,  # Statement 3
    lambda t, f: f == 3,  # Statement 4
    lambda t, f: (t == 4) != (f == 3),  # Statement 5
    lambda t, f: t in {2, 3, 5, 7},  # Statement 6
    lambda t, f: f in {4, 6},  # Statement 7
]

# Check all combinations of true/false statements
valid_solutions = 0
for true_statements in range(8):  # 0 to 7 true statements
    false_statements = 7 - true_statements
    if all(statement(true_statements, false_statements) for statement in statements):
        valid_solutions += 1

print(valid_solutions)