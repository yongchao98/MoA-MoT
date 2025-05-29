# Define the truth values for R, Q, and P
truth_values = [(False, False, False), (False, False, True), (False, True, False), (False, True, True),
                (True, False, False), (True, False, True), (True, True, False), (True, True, True)]

# Check which combinations satisfy the premises
valid_combinations = []
for R, Q, P in truth_values:
    premise1 = not (R and Q)
    premise2 = (Q == P)
    if premise1 and premise2:
        valid_combinations.append((R, Q, P))

# Determine the conclusion
conclusion = all(R == False or Q == False for R, Q, P in valid_combinations)

print(conclusion)