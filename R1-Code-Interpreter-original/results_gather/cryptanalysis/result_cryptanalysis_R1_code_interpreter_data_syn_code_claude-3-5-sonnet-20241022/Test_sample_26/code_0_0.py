from itertools import permutations

def check_guess(guess, actual):
    # Check if a guess matches all the feedback conditions
    numbers_match = [False, False]
    letters_match = [False, False]
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            numbers_match[i] = True
    
    # Check letters
    for i in range(2, 4):
        if guess[i] == actual[i]:
            letters_match[i-2] = True
            
    return numbers_match, letters_match

def is_valid_solution(solution):
    # Convert solution to string format for easier comparison
    sol = [str(solution[0]), str(solution[1]), solution[2], solution[3]]
    
    # Check all conditions
    
    # Condition 1: 83VM
    nm, lm = check_guess(['8','3','V','M'], sol)
    if any(nm) or (lm != [False, False]) or (sol[2:].count('V') + sol[2:].count('M') != 1):
        return False
    
    # Condition 2: 40WL
    if 'W' in sol[2:] or 'L' in sol[2:]:
        return False
    
    # Condition 3: 50VD
    if 'D' in sol[2:]:
        return False
    
    # Condition 4: 40MC
    if 'M' in sol[2:] or 'C' in sol[2:]:
        return False
    
    # Condition 5: 18JC
    nm, lm = check_guess(['1','8','J','C'], sol)
    if sum(nm) != 1 or any(lm):
        return False
    if sol[0] != '1' and int(sol[0]) > 1:  # First number must be 1 or less than 1
        return False
        
    # Condition 6: 41TW
    nm, lm = check_guess(['4','1','T','W'], sol)
    if sum(nm) != 1 or any(lm):
        return False
    
    # Condition 7: 07BF
    if 'B' in sol[2:]:
        return False
    
    return True

# Generate all possible combinations
numbers = range(10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

valid_solutions = []

# Try all possible combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        solution = [n1, n2, l1, l2]
                        if is_valid_solution(solution):
                            valid_solutions.append(solution)

print(valid_solutions)