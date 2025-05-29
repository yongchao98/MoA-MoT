from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K, L, M = perm
    
    # Check all conditions
    if B + L != 35:  # Condition 1
        return False
    if J <= D:  # Condition 2
        return False
    if L - J != -8:  # Condition 3
        return False
    if abs(I - 3.6 * H) > 0.1:  # Condition 4 (using small tolerance)
        return False
    if B + I != 64:  # Condition 5
        return False
    if H + I != 46:  # Condition 6
        return False
    if C + L != 31:  # Condition 7
        return False
    if G + K != 19:  # Condition 8
        return False
    if M - J != 30:  # Condition 9
        return False
    if abs(K - 3.2 * E) > 0.1:  # Condition 10 (using small tolerance)
        return False
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break