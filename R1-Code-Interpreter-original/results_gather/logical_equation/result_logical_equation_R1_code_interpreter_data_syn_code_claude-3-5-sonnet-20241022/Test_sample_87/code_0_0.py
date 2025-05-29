from itertools import permutations

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

def check_conditions(values):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if not all([
        I - J == -21,           # 1
        J - D == 14,            # 2
        K == 3.5 * H,          # 3
        abs(E - 3.2 * C) < 0.1,  # 4 (using tolerance for floating point)
        J - H == 22,            # 5
        M - I == 33,            # 6
        abs(J - 4.8 * C) < 0.1,  # 7
        abs(E - 1.6 * D) < 0.1,  # 8
        I - L == -12,           # 9
        abs(M - 2.4 * L) < 0.1,  # 10
        C - E == -11            # 11
    ]):
        return False
    return True

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break