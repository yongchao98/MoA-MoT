from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if D - K != 23: return False
    if L != 3.6 * A: return False
    if A + H != 25: return False
    if K != 2.5 * C: return False
    if B <= F: return False
    if E != 2.4 * A: return False
    if H != 1.5 * A: return False
    if I - G != -42: return False
    if I - E != -21: return False
    if L != 1.5 * E: return False
    
    return True

# Known values from our analysis
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
A = 10
H = 15
E = 24
I = 3
G = 45
L = 36

# Find remaining values
remaining_nums = [n for n in numbers if n not in [A, H, E, I, G, L]]
remaining_positions = [1, 2, 3, 4, 5, 11, 12]  # Indices for B, C, D, F, J, K, M

for perm in permutations(remaining_nums):
    values = [0] * 13
    values[0] = A
    values[7] = H
    values[4] = E
    values[8] = I
    values[6] = G
    values[11] = L
    
    # Fill in the remaining values
    for pos, val in zip(remaining_positions, perm):
        values[pos] = val
    
    if check_conditions(values):
        print(values)
        break