from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if E + K != 81: return False
    if K <= J: return False
    if abs(E - 1.5*G) > 0.0001: return False
    if B + J != 5: return False
    if M - F != -6: return False
    if abs(K - 4.5*M) > 0.0001: return False
    if M <= B: return False
    if abs(F - 3.2*C) > 0.0001: return False
    if C - E != -31: return False
    if abs(L - 1.4*C) > 0.0001: return False
    if A + M != 38: return False
    if F + I != 112: return False
    
    return True

numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all possible permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        break