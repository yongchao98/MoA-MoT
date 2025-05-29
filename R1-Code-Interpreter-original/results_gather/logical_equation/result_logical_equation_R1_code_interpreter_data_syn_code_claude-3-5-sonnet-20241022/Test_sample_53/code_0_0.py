from itertools import permutations

# The available numbers
numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

def check_constraints(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all constraints
    if F - K != 36: return False
    if B - K != 111: return False
    if B + E != 230: return False
    if H + J != 25: return False
    if E - A != 77: return False
    if H + I != 36: return False
    if C + E != 305: return False
    if C - I != 205: return False
    if C + J != 234: return False
    if G != 3.0 * H: return False
    if A + K != 42: return False
    
    return True

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm):
        # Print in required format
        print(f"<<<{list(perm)}>>>")
        break