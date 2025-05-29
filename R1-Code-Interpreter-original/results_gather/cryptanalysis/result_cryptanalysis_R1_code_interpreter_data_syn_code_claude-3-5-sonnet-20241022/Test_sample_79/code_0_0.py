from itertools import permutations

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    
    # Condition from guess 3: one number from 7,8 is correct but wrong position
    # and one is too large
    if not ((('7' in numbers and numbers.index('7') != 0) or 
             ('8' in numbers and numbers.index('8') != 1)) and 
            not ('7' in numbers and '8' in numbers)):
        return False
    
    # All these number combinations must be wrong
    wrong_numbers = ['85', '19', '78', '61', '41', '12']
    if ''.join(numbers) in wrong_numbers:
        return False
    
    # Letters conditions
    if 'U' in letters or 'J' in letters:  # from guess 1
        return False
    
    # From guess 2: B or P correct but wrong position
    if not (('B' in letters and letters.index('B') != letters.index('B')) or 
            ('P' in letters and letters.index('P') != 1)):
        return False
        
    # From guess 3: L or U correct and in position
    if not ('L' in letters):  # We know U is wrong from earlier
        return False
        
    # V, Z, P, R are too late in alphabet
    if any(c in letters for c in 'VZPR'):
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

valid_solutions = []
for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        guess = [n1, n2, l1, l2]
        if check_conditions(guess):
            valid_solutions.append(guess)

print(valid_solutions)