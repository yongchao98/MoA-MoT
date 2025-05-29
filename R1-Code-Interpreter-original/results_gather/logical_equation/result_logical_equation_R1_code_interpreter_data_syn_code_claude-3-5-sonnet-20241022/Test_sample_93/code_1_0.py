def is_valid(assignment, pos, val, values):
    # Convert assignment to list for easier indexing
    current = [-1] * 13
    for k, v in assignment.items():
        current[k] = v
    current[pos] = val
    
    # Only check constraints where we have all necessary values assigned
    if current[2] != -1 and current[4] != -1:  # C + E = 141
        if current[2] + current[4] != 141: return False
    if current[9] != -1 and current[11] != -1:  # J + L = 33
        if current[9] + current[11] != 33: return False
    if current[10] != -1 and current[7] != -1:  # K = 2.4H
        if abs(current[10] - 2.4*current[7]) > 0.1: return False
    if current[1] != -1 and current[9] != -1:  # B - J = 10
        if current[1] - current[9] != 10: return False
    if current[4] != -1 and current[2] != -1:  # E - C = 51
        if current[4] - current[2] != 51: return False
    if current[2] != -1 and current[5] != -1:  # C > F
        if current[2] <= current[5]: return False
    if current[12] != -1 and current[9] != -1:  # M = 3.2J
        if abs(current[12] - 3.2*current[9]) > 0.1: return False
    if current[4] != -1 and current[10] != -1:  # E + K = 120
        if current[4] + current[10] != 120: return False
    if current[9] != -1 and current[7] != -1:  # J - H = -5
        if current[9] - current[7] != -5: return False
    if current[10] != -1 and current[1] != -1:  # K - B = 9
        if current[10] - current[1] != 9: return False
    if current[1] != -1 and current[12] != -1:  # B - M = -1
        if current[1] - current[12] != -1: return False
    if current[2] != -1 and current[12] != -1:  # C > M
        if current[2] <= current[12]: return False
    return True

def solve():
    values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    assignment = {}
    used = set()
    
    # Start with most constrained variables
    order = [2, 4, 9, 7, 10, 1, 12, 5, 11, 0, 3, 6, 8]  # C,E,J,H,K,B,M,F,L,A,D,G,I
    
    def backtrack(pos):
        if pos == 13:
            return True
        
        idx = order[pos]
        for val in values:
            if val in used:
                continue
            if is_valid(assignment, idx, val, values):
                assignment[idx] = val
                used.add(val)
                if backtrack(pos + 1):
                    return True
                used.remove(val)
                del assignment[idx]
        return False
    
    if backtrack(0):
        result = [-1] * 13
        for k, v in assignment.items():
            result[k] = v
        print(result)
    else:
        print("No solution found")

solve()