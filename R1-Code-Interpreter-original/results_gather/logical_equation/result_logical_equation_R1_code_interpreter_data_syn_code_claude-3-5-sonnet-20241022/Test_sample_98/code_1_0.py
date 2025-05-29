def check_solution(numbers):
    # Pre-assign known values
    solution = [0] * 13  # A to M
    solution[9] = 36  # J = 36
    solution[10] = 45  # K = 45
    solution[0] = 96  # A = 96
    solution[1] = 50  # B = 50
    
    remaining_numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28]
    remaining_positions = [2, 3, 4, 5, 6, 7, 8, 11, 12]
    
    def is_valid(temp_solution):
        A, B, C, D, E, F, G, H, I, J, K, L, M = temp_solution
        
        # Check all constraints
        if A + F != 30:  # A = 96, so F = 16
            return False
        if C + D != 31:
            return False
        if G - L != 45:
            return False
        if L - C != -19:
            return False
        if abs(I - 1.5*K) > 0.1:  # I should be close to 67.5
            return False
        if A <= D:
            return False
        if abs(C - 1.6*I) > 0.1:
            return False
        if B + G != 86:  # B = 50, so G = 36
            return False
        
        return True
    
    from itertools import permutations
    for perm in permutations(remaining_numbers):
        temp_solution = solution.copy()
        for pos, num in zip(remaining_positions, perm):
            temp_solution[pos] = num
        if is_valid(temp_solution):
            return temp_solution
    
    return None

result = check_solution([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96])
if result:
    print(result)