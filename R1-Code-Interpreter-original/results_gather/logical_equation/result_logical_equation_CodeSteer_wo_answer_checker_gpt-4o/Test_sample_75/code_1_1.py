from itertools import permutations

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Define a function to check if a given assignment satisfies all constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    # Check all constraints
    if not (L - I == -23): return False
    if not (I + M == 44): return False
    if not (D + K == 98): return False
    if not (F + I == 73): return False
    if not (F > J): return False
    if not (B - C == 26): return False
    if not (B == 2.4 * A): return False
    if not (C + F == 55): return False
    if not (F - J == 38): return False
    if not (H - F == 5): return False
    if not (J == 1.4 * L): return False
    
    return True

# Try all permutations of the possible values
for perm in permutations(possible_values):
    if satisfies_constraints(perm):
        print(f"<<<{list(perm)}>>>")
        break