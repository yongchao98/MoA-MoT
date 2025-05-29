from itertools import permutations

def check_equations(values, letters):
    A, B, C, D, E, F, G, H, I = values
    
    # Check all equations
    if D != 1.5 * I:  # D = 1.5I
        return False
    if I != 2.5 * H:  # I = 2.5H
    	return False
    if C + H != 52:   # C + H = 52
        return False
    if I != 4.0 * C:  # I = 4.0C
        return False
    if C - H != -12:  # C - H = -12
        return False
    if D + I != 200:  # D + I = 200
        return False
    if B - H != -23:  # B - H = -23
        return False
    if H != 1.6 * C:  # H = 1.6C
        return False
    
    return True

# Numbers to assign
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

# Try all permutations
for perm in permutations(numbers):
    if check_equations(perm, letters):
        # Create a dictionary for easy lookup
        solution = dict(zip(letters, perm))
        # Print in alphabetical order
        result = [solution[letter] for letter in sorted(solution.keys())]
        print(result)
        break