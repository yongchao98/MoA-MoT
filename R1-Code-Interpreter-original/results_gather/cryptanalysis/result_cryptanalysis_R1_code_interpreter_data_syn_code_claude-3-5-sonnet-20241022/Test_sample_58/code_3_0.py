def check_final_solution(guess):
    # Guess format: [num1, num2, letter1, 'Z']
    
    # First, check if second letter is Z (we know this is correct)
    if guess[3] != 'Z':
        return False
    
    # First letter can't be Z (no repeating letters allowed)
    if guess[2] == 'Z':
        return False
    
    # Check number positions and values
    # From guess #2 (79EN): 7 is correct but wrong position
    # From guess #3 (70GE): 7 is correct but wrong position
    # This means 7 must be in first position
    if guess[0] != '7':
        return False
    
    # From guess #5 (05LE): 5 is correct but wrong position
    # Combined with above, 5 must be in second position
    if guess[1] != '5':
        return False
    
    # Letter conditions
    letter_val = ord(guess[2])
    
    # Must be after T due to multiple "too early" feedbacks
    if letter_val <= ord('T'):
        return False
    
    # Must be Y because:
    # 1. We know it's after T
    # 2. From guess #7, one letter is Z (which we know is second position)
    # 3. The other letter in guess #7 (F) is too early
    # 4. Therefore, our first letter must be Y to satisfy all conditions
    if guess[2] != 'Y':
        return False
    
    return True

candidates = [['7', '5', l1, 'Z'] for l1 in 'UVWXYZ'] + [['5', '7', l1, 'Z'] for l1 in 'UVWXYZ']
valid_solutions = [cand for cand in candidates if check_final_solution(cand)]
print(valid_solutions)