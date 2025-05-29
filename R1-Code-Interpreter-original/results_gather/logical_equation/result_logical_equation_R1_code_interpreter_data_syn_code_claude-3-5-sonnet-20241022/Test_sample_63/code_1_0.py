def verify_and_print(solution):
    valid_numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
    A, B, C, D, E, F, G, H, I, J, K, L, M = solution
    
    # Check if all numbers are from valid set
    if not all(x in valid_numbers for x in solution):
        return False
        
    # Check all conditions
    if A - E != 4: return False
    if E - A != -4: return False
    if I != 4.0 * J: return False
    if D - B != -14: return False
    if H + M != 26: return False
    if abs(L - 3.2 * A) > 0.1: return False
    if abs(L - 1.6 * F) > 0.1: return False
    if F + L != 26: return False
    if H - I != -4: return False
    if abs(D - 1.5 * H) > 0.1: return False
    
    print(solution)
    return True

# Let's solve with known constraints:
H = 16  # Start with H as it's in multiple equations
D = 24  # D = 1.5H
I = 28  # H - I = -4
J = 7   # I = 4J
B = 36  # D - B = -14
M = 10  # H + M = 26
F = 10  # Try F = 10
L = 16  # F + L = 26
A = 5   # L ≈ 3.2A
E = 1   # A - E = 4
# Remaining numbers: [2, 3, 15, 45, 50]
C = 3
G = 2
K = 15

solution = [A, B, C, D, E, F, G, H, I, J, K, L, M]
verify_and_print(solution)