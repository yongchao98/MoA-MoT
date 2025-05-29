from itertools import permutations

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if D - A != 43: return False
    if abs(F - 1.5 * C) > 0.01: return False  # Using small epsilon for float comparison
    if A + F != 43: return False
    if abs(C - 1.6 * L) > 0.01: return False
    if abs(L - 3.0 * J) > 0.01: return False
    if F - J != 31: return False
    if H <= L: return False
    if abs(M - 4.5 * G) > 0.01: return False
    if E + G != 12: return False
    if J + K != 33: return False
    if abs(F - 3.6 * G) > 0.01: return False
    if M <= A: return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Try all permutations
for p in permutations(numbers):
    if check_conditions(p):
        print(list(p))
        break