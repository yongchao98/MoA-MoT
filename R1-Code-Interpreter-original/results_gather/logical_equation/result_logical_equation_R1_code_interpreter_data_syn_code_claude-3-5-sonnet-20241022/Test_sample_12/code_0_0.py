from itertools import permutations

def check_conditions(values, letters):
    # Create a dictionary for easy reference
    d = dict(zip(letters, values))
    
    # Check all conditions
    if d['D'] + d['H'] != 36:  # Condition 1
        return False
    if d['G'] - d['I'] != -79:  # Condition 2
        return False
    if d['C'] - d['B'] != -20:  # Condition 3
        return False
    if d['D'] != 3.0 * d['E']:  # Condition 4
        return False
    if d['B'] != 4.0 * d['D']:  # Condition 5
        return False
    if d['I'] <= d['C']:  # Condition 6
        return False
    if d['B'] + d['C'] != 52:  # Condition 7
        return False
    return True

# Available numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']

# Try all permutations
for perm in permutations(numbers):
    if check_conditions(perm, letters):
        print(list(perm))
        break