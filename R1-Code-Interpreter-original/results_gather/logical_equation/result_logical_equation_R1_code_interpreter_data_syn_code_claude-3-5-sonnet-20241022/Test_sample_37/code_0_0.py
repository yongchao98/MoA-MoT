from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K = [values[letters.index(x)] for x in 'ABCDEFGHIJK']
    
    # Check all constraints
    if I != 1.5 * H:  # 1
        return False
    if F != 4.0 * B:  # 2
        return False
    if C + E != 169:  # 3
        return False
    if E != 3.0 * J:  # 4
        return False
    if J - G != -11:  # 5
        return False
    if C <= F:  # 6
        return False
    if I + K != 305:  # 7
        return False
    if I - D != 222:  # 8
        return False
    if J - F != -20:  # 9
        return False
    if F + H != 186:  # 10
        return False
    
    return True

def find_solution():
    numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
    letters = list('ABCDEFGHIJK')
    
    for perm in permutations(numbers):
        if check_constraints(perm, letters):
            return list(perm)
    
    return None

solution = find_solution()
print(solution)