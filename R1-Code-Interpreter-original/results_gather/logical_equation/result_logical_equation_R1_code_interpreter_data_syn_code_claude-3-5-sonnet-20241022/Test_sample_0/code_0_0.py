def check_conditions(perm):
    A, B, C, D, E, F, G, H, I = perm
    
    # Check all conditions
    if A - C != -44: return False
    if C <= H: return False
    if I != 3 * D: return False
    if B - C != -79: return False
    if D != 3 * B: return False
    if A != 4 * I: return False
    if B - E != -3: return False
    
    return True

def find_solution():
    numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
    
    from itertools import permutations
    for perm in permutations(numbers):
        if check_conditions(perm):
            return list(perm)
    
    return None

solution = find_solution()
print(solution)