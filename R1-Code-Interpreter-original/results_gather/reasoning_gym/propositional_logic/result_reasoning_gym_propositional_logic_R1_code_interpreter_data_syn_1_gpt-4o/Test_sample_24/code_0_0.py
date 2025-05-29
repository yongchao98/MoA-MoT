from itertools import product

# Define the premises
def premise1(R, S):
    return R == S

def premise2(R, P):
    return (R and P) or not R

# Define the conclusion
def conclusion(R, S, P):
    return not R or (S and P)

# Check all combinations of truth values for R, S, P
truth_values = [True, False]
valid_conclusion = True

for R, S, P in product(truth_values, repeat=3):
    if premise1(R, S) and premise2(R, P):
        if not conclusion(R, S, P):
            valid_conclusion = False
            break

print(valid_conclusion)