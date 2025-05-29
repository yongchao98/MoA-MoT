from itertools import permutations

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if A + I != 300:  # Condition 1
        return False
    if abs(E - 2.4*C) > 0.1:  # Condition 2 (using small tolerance)
        return False
    if abs(G - 3.0*K) > 0.1:  # Condition 3
        return False
    if J <= A:  # Condition 4
        return False
    if abs(A - G - 66) > 0.1:  # Condition 5
        return False
    if abs(J - 4.0*C) > 0.1:  # Condition 6
        return False
    if abs(I - 3.0*A) > 0.1:  # Condition 7
        return False
    if abs(F - 2.0*A) > 0.1:  # Condition 8
        return False
    if abs(I - 1.5*F) > 0.1:  # Condition 9
        return False
    if abs(A - C - 55) > 0.1:  # Condition 10
        return False
    if abs(H - A + 59) > 0.1:  # Condition 11
        return False
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm):
        print(list(perm))
        break