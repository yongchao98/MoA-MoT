from itertools import permutations

def check_conditions(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all conditions
    if F + H != 6:  # Condition 1
        return False
    if abs(E - 2.4*D) > 0.01:  # Condition 2
        return False
    if abs(B - 2.0*F) > 0.01:  # Condition 3
        return False
    if abs(C - 3.0*M) > 0.01:  # Condition 4
        return False
    if I - B != 34:  # Condition 5
        return False
    if abs(H - 2.5*B) > 0.01:  # Condition 6
        return False
    if E + M != 39:  # Condition 7
        return False
    if abs(L - 4.0*J) > 0.01:  # Condition 8
        return False
    if abs(M - 3.0*H) > 0.01:  # Condition 9
        return False
    if abs(L - 2.8*D) > 0.01:  # Condition 10
        return False
    
    return True

numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break