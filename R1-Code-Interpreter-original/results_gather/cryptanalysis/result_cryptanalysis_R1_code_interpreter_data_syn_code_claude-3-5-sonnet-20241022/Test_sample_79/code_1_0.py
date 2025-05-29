from itertools import permutations

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    
    # Check number conditions
    # From guess 3: one of 7,8 correct but wrong position, other too large
    has_valid_7_8 = False
    if '7' in numbers and '8' not in numbers:
        if numbers.index('7') == 1:  # 7 in wrong position
            has_valid_7_8 = True
    elif '8' in numbers and '7' not in numbers:
        if numbers.index('8') == 0:  # 8 in wrong position
            has_valid_7_8 = True
    if not has_valid_7_8:
        return False
    
    # All these combinations must be wrong
    wrong_numbers = ['85', '19', '78', '61', '41', '12']
    if ''.join(numbers) in wrong_numbers:
        return False
    
    # Letter conditions
    # From guess 1: U and J are wrong
    if 'U' in letters or 'J' in letters:
        return False
    
    # From guess 2: B or P correct but wrong position
    has_b_p_wrong_pos = False
    if 'B' in letters:
        if letters.index('B') != 1:  # B not in same position as in guess 2
            has_b_p_wrong_pos = True
    elif 'P' in letters:
        if letters.index('P') != 1:
            has_b_p_wrong_pos = True
    if not has_b_p_wrong_pos:
        return False
    
    # From guess 3: L or U correct and in correct position
    # We know U is wrong, so L must be correct and in position
    if 'L' not in letters:
        return False
    
    # V, Z, P, R are too late in alphabet
    if any(c in letters for c in 'VZPR'):
        return False
    
    # From guess 5: E or L correct but wrong position
    if 'E' in letters:
        if letters.index('E') == 0:  # E must be in wrong position
            return False
    elif 'L' in letters:
        if letters.index('L') == 0:  # L must be in wrong position
            return False
    
    return True

# Generate all possible combinations
valid_solutions = []
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 76)]  # A-L only, as we know later letters are too late

for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        guess = list(nums + lets)
        if check_conditions(guess):
            valid_solutions.append(guess)

print(valid_solutions)