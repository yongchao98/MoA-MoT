from itertools import permutations

# Available numbers
numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]

def check_conditions(perm):
    A, B, C, D, E, F, G, H, I, J, K = perm
    
    # Check all conditions
    if F - C != -141:  # Condition 1
        return False
    if C <= D:  # Condition 2
        return False
    if E + K != 305:  # Condition 3
        return False
    if abs(K - 3.0*A) > 0.1:  # Condition 4
        return False
    if abs(H - 3.0*F) > 0.1:  # Condition 5
        return False
    if I + J != 271:  # Condition 6
        return False
    if abs(K - 1.5*C) > 0.1:  # Condition 7
        return False
    if K - A != 150:  # Condition 8
        return False
    if F + G != 123:  # Condition 9
        return False
    if B + C != 153:  # Condition 10
        return False
    
    return True

# Try all permutations and print in required format
for p in permutations(numbers):
    if check_conditions(p):
        print(f"<<<{list(p)}>>>")
        break