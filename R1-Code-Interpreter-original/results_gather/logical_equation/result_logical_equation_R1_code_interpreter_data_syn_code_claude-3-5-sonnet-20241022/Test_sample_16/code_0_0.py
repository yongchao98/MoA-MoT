from itertools import permutations

def check_equations(values, letters):
    # Create a dictionary mapping letters to values
    d = dict(zip(letters, values))
    
    # Check all equations
    if d['C'] != 3.0 * d['D']:  # C = 3.0D
        return False
    if d['C'] != 1.5 * d['F']:  # C = 1.5F
        return False
    if d['D'] - d['C'] != -6:   # D - C = -6
        return False
    if d['F'] != 2.0 * d['D']:  # F = 2.0D
        return False
    if d['D'] + d['E'] != 124:  # D + E = 124
        return False
    if d['A'] + d['H'] != 112:  # A + H = 112
        return False
    if d['I'] != 1.5 * d['A']:  # I = 1.5A
        return False
    if d['D'] + d['G'] != 23:   # D + G = 23
        return False
    return True

# Available numbers
numbers = [3, 6, 9, 20, 32, 36, 80, 120, 121]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

# Try all permutations
for perm in permutations(numbers):
    if check_equations(perm, letters):
        # Found the solution
        result = list(perm)
        print(result)
        break