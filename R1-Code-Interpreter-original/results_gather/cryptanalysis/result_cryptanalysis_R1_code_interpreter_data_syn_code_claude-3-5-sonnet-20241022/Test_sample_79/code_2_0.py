from itertools import permutations

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    n1, n2 = map(int, numbers)
    l1, l2 = letters
    
    # Numbers that can't be used (from all guesses)
    forbidden_numbers = {1, 2, 4, 5, 6, 8, 9}
    if n1 in forbidden_numbers or n2 in forbidden_numbers:
        return False
    
    # From guess 3: one of 7,8 correct but wrong position
    # Since 8 is forbidden, 7 must be used but in wrong position
    if n1 == 7:  # 7 must be in second position
        return False
    if n2 != 7:  # 7 must be present
        return False
    
    # Letters that are too late in alphabet
    if any(c in letters for c in 'JPRUVZ'):
        return False
    
    # From guess 2: B or P correct but wrong position
    # Since P is too late, B must be present but in different position than in guess 2
    if 'B' not in letters:
        return False
    if letters.index('B') == 1:  # B can't be in same position as in guess 2
        return False
    
    # From guess 3: L must be in correct position (position 3)
    if letters[1] != 'L':  # L must be in second position
        return False
    
    # From guess 5: E or L correct but wrong position
    # We know L is correct in position, so E must not be used
    if 'E' in letters:
        return False
    
    return True

valid_solutions = []
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 77)]  # A-L only

for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        guess = [n1, n2, l1, l2]
        if check_conditions(guess):
            valid_solutions.append(guess)

print(valid_solutions)