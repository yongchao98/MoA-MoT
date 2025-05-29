from itertools import permutations

def find_solution():
    numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
    
    # Pre-filter possible values based on direct relationships
    for J in numbers:
        # From equation 1: I = J - 21
        I = J - 21
        if I not in numbers:
            continue
            
        # From equation 2: D = J - 14
        D = J - 14
        if D not in numbers:
            continue
            
        # From equation 5: H = J - 22
        H = J - 22
        if H not in numbers:
            continue
            
        # From equation 9: L = I + 12
        L = I + 12
        if L not in numbers:
            continue
            
        # From equation 3: K = 3.5H
        K = int(3.5 * H)
        if K not in numbers:
            continue
            
        # From equation 6: M = I + 33
        M = I + 33
        if M not in numbers:
            continue
            
        # From equation 7: C = J/4.8
        C = round(J/4.8, 1)
        if int(C) not in numbers:
            continue
            
        # From equation 8 and 4: E = 1.6D and E = 3.2C
        E = round(1.6 * D, 1)
        if int(E) not in numbers:
            continue
            
        if abs(E - 3.2 * C) > 0.1:
            continue
            
        # Check if C - E = -11
        if C - E != -11:
            continue
            
        # Find remaining values (A, B, F, G)
        remaining = set(numbers) - {int(x) for x in [C, D, E, H, I, J, K, L, M]}
        
        for A, B, F, G in permutations(remaining):
            result = [A, B, int(C), int(D), int(E), F, G, int(H), int(I), int(J), int(K), int(L), int(M)]
            print(result)
            return

find_solution()