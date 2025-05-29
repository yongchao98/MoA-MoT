from itertools import permutations

def check_equations(values, letters):
    A, B, C, D, E, F, G, H, I, J, K = values
    
    # Check all equations
    if E != 1.5 * K: return False
    if D + I != 155: return False
    if K != 2.0 * D: return False
    if E != 3.0 * D: return False
    if D - C != -46: return False
    if F + H != 51: return False
    if H + I != 83: return False
    if F != 3.0 * G: return False
    if E + K != 375: return False
    if F - D != -27: return False
    
    return True

numbers = [3, 9, 16, 27, 36, 48, 75, 80, 121, 150, 225]
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Try all permutations
for perm in permutations(numbers):
    if check_equations(perm, letters):
        # Create a dictionary mapping letters to values
        result = dict(zip(letters, perm))
        # Sort by letters and get values
        sorted_result = [result[letter] for letter in sorted(result.keys())]
        print(sorted_result)
        break