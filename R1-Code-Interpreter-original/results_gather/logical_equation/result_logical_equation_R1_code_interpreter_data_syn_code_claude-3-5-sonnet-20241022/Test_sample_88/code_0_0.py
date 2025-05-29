from itertools import permutations

def check_constraints(values, letters):
    A, B, C, D, E, F, G, H, I, J, K, L, M = values
    
    # Check all constraints
    if J != 3.5 * D:  # J = 3.5D
        return False
    if K - F != -12:  # K - F = -12
        return False
    if G + M != 26:  # G + M = 26
        return False
    if J + K != 31:  # J + K = 31
        return False
    if C - D != 13:  # C - D = 13
        return False
    if G != 1.6 * M:  # G = 1.6M
        return False
    if L <= G:  # L > G
        return False
    if A + L != 124:  # A + L = 124
        return False
    if M - A != -18:  # M - A = -18
        return False
    if C != 3.0 * H:  # C = 3.0H
        return False
    if J <= D:  # J > D
        return False
    
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_constraints(perm, letters):
        # Create a dictionary for easy lookup
        solution = dict(zip(letters, perm))
        # Print in alphabetical order
        result = [solution[letter] for letter in sorted(letters)]
        print(result)
        break