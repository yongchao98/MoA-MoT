from itertools import permutations

def is_valid_partial(B, D, E, F, H, I, K, L, M):
    if B is not None and I is not None and B + I != 31: return False
    if I is not None and D is not None and I != 1.5 * D: return False
    if K is not None and E is not None and K - E != -2: return False
    if B is not None and H is not None and B + H != 57: return False
    if F is not None and I is not None and F != 4.0 * I: return False
    if K is not None and D is not None and K - D != -13: return False
    if M is not None and E is not None and M != 2.0 * E: return False
    if I is not None and L is not None and I + L != 69: return False
    if M is not None and B is not None and M <= B: return False
    if D is not None and H is not None and D - H != -34: return False
    return True

# First, let's determine some relationships
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# We know I = 1.5D, so I must be 1.5 times some number in the list
possible_I = []
possible_D = []
for n in numbers:
    if n * 1.5 in numbers:
        possible_D.append(n)
        possible_I.append(n * 1.5)

# We know F = 4.0I
possible_F = []
for i in possible_I:
    if i * 4.0 in numbers:
        possible_F.append(i * 4.0)

# We know M = 2.0E
possible_M = []
possible_E = []
for n in numbers:
    if n * 2.0 in numbers:
        possible_E.append(n)
        possible_M.append(n * 2.0)

result = None
for D in possible_D:
    I = D * 1.5
    F = I * 4.0
    if F not in numbers:
        continue
    
    for E in possible_E:
        M = E * 2.0
        K = D - 13
        
        if K not in numbers or M not in numbers:
            continue
            
        remaining_nums = set(numbers) - {D, I, F, E, M, K}
        
        for B in remaining_nums:
            if B + I != 31:
                continue
            if B + H != 57:
                continue
                
            H = 57 - B
            if H not in remaining_nums:
                continue
                
            if D - H != -34:
                continue
                
            remaining_nums2 = remaining_nums - {B, H}
            
            for L in remaining_nums2:
                if I + L != 69:
                    continue
                    
                remaining_nums3 = remaining_nums2 - {L}
                
                for A, G, J in permutations(remaining_nums3, 3):
                    if A + G != 43:
                        continue
                        
                    result = [A, B, C, D, E, F, G, H, I, J, K, L, M]
                    break
                if result:
                    break
            if result:
                break
        if result:
            break
    if result:
        break

if result:
    print(result)