from itertools import permutations

# Given numbers
numbers = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Constraints as functions
def constraints(A, B, C, D, E, F, G, H, I, J, K, L, M):
    return (M == 2.4 * L and
            F == 3.0 * L and
            E + I == 8 and
            K - H == -8 and
            H > I and
            L + M == 51 and
            C + J == 9 and
            G == 2.0 * I and
            A - D == 49 and
            H == 1.5 * K)

# Backtracking function
def backtrack(assignment, index):
    if index == len(letters):
        # If all letters are assigned, check constraints
        if constraints(*assignment):
            print(f"<<<{assignment}>>>")
            return True
        return False
    
    for number in numbers:
        if number not in assignment:
            assignment[index] = number
            if backtrack(assignment, index + 1):
                return True
            assignment[index] = None
    return False

# Letters to assign
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Initialize assignment
assignment = [None] * len(letters)

# Start backtracking
backtrack(assignment, 0)