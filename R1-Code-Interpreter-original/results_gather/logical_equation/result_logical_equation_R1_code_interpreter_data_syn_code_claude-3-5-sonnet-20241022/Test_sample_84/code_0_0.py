from itertools import permutations

def check_conditions(values, letters):
    # Create a dictionary for easy lookup
    d = dict(zip(letters, values))
    
    # Check all conditions
    if abs(d['B'] - 1.4 * d['G']) > 0.1: return False  # B = 1.4G
    if abs(d['A'] - 3.2 * d['G']) > 0.1: return False  # A = 3.2G
    if abs(d['J'] - 3.0 * d['G']) > 0.1: return False  # J = 3.0G
    if abs(d['A'] - 1.6 * d['K']) > 0.1: return False  # A = 1.6K
    if abs(d['A'] + d['J'] - 31) > 0.1: return False   # A + J = 31
    if abs(d['A'] + d['B'] - 23) > 0.1: return False   # A + B = 23
    if abs(d['F'] - 2.4 * d['K']) > 0.1: return False  # F = 2.4K
    if abs(d['M'] - d['B'] - 89) > 0.1: return False   # M - B = 89
    if abs(d['I'] - d['B'] - 43) > 0.1: return False   # I - B = 43
    if abs(d['F'] - 4.8 * d['G']) > 0.1: return False  # F = 4.8G
    if abs(d['B'] - 3.5 * d['L']) > 0.1: return False  # B = 3.5L
    return True

# Available numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        result = list(perm)
        print(result)
        break