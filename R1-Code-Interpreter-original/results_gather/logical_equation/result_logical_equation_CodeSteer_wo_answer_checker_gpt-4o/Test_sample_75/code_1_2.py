from itertools import permutations

# Define the possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Pre-calculate some values based on the equations
# From L - I = -23, we have L = I - 23
# From I + M = 44, we have M = 44 - I
# From F + I = 73, we have F = 73 - I
# From F - J = 38, we have J = F - 38
# From H - F = 5, we have H = F + 5

def satisfies_constraints(A, B, C, D, E, F, G, H, I, J, K, L, M):
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

# Reduce the number of permutations by solving some equations first
for I in possible_values:
    L = I - 23
    M = 44 - I
    F = 73 - I
    J = F - 38
    H = F + 5
    
    # Check if L, M, F, J, H are in possible values
    if L in possible_values and M in possible_values and F in possible_values and J in possible_values and H in possible_values:
        # Remove used values from possible values
        remaining_values = [v for v in possible_values if v not in {I, L, M, F, J, H}]
        
        # Try permutations of the remaining values for A, B, C, D, E, G, K
        for perm in permutations(remaining_values):
            A, B, C, D, E, G, K = perm
            if satisfies_constraints(A, B, C, D, E, F, G, H, I, J, K, L, M):
                print(f"<<<{[A, B, C, D, E, F, G, H, I, J, K, L, M]}>>>")
                break