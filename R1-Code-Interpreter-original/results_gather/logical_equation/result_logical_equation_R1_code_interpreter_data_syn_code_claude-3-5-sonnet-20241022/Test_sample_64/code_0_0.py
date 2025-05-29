from itertools import permutations

def check_constraints(values, nums):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if E + H != 31: return False
    if F - K != 4: return False
    if D + K != 31: return False
    if abs(I - 3.6 * L) > 0.0001: return False  # Using small epsilon for float comparison
    if abs(J - 1.6 * E) > 0.0001: return False
    if C - I != -35: return False
    if A - H != -14: return False
    if abs(I - 2.4 * E) > 0.0001: return False
    if D + G != 33: return False
    if J <= G: return False
    
    return True

# The available numbers
nums = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all permutations
for p in permutations(nums):
    if check_constraints(p, nums):
        print(list(p))
        break