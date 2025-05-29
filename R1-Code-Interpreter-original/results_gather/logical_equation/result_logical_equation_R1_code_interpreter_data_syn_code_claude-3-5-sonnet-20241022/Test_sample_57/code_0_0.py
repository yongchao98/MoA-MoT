from itertools import permutations

def check_constraints(values, positions):
    A, B, C, D, E, F, G, H, I, J, K = [values[positions[x]] for x in 'ABCDEFGHIJK']
    
    # Check all equality constraints
    if C - G != -36: return False
    if C - F != -17: return False
    if B + F != 68: return False
    if E + H != 137: return False
    if A + K != 375: return False
    if C + H != 124: return False
    if F - G != -19: return False
    if C + G != 42: return False
    
    # Check inequality constraints
    if not (A > F): return False
    if not (J > C): return False
    if not (A > C): return False
    
    return True

numbers = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]
letters = 'ABCDEFGHIJK'

# Try all possible permutations
for perm in permutations(range(11)):
    positions = {letter: pos for letter, pos in zip(letters, perm)}
    if check_constraints(numbers, positions):
        result = [numbers[positions[letter]] for letter in letters]
        print(result)
        break