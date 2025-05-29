from itertools import permutations

def check_equations(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all equations
    if not abs(D - 3.0 * B) < 0.1:  # D = 3.0B
        return False
    if A + F != 52:  # A + F = 52
        return False
    if D - H != -71:  # D - H = -71
        return False
    if B + I != 30:  # B + I = 30
        return False
    if D + E != 13:  # D + E = 13
        return False
    if not abs(I - 3.0 * D) < 0.1:  # I = 3.0D
        return False
    if not abs(F - 4.0 * E) < 0.1:  # F = 4.0E
        return False
    
    return True

numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
found = False

for perm in permutations(numbers):
    if check_equations(perm):
        print(list(perm))
        found = True
        break

if not found:
    print("No solution found")