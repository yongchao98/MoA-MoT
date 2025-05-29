from itertools import product

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    
    # Check number conditions
    if int(numbers) <= 2: return False  # too small
    if int(numbers) >= 97: return False  # too large
    
    # Known incorrect numbers
    wrong_numbers = ['97','70','17','72','02','54','40','09']
    if numbers in wrong_numbers: return False
    
    # Check 96CL condition - one number must be 6 or 9 in correct position
    if '6' not in numbers[1] and '9' not in numbers[0]: return False
    
    # Check letter conditions
    wrong_letters = 'YSGFXPCLEW'
    if letters[0] in wrong_letters or letters[1] in wrong_letters: return False
    
    # Letters can't be too early (before G) or too late (after R)
    for letter in letters:
        if letter < 'G' or letter > 'R': return False
    
    # From guess 11, X or R must be in correct position
    if letters[0] != 'X' and letters[1] != 'R': return False
    
    return True

# Generate all possible combinations
digits = '0123456789'
letters = 'GHIJKLMNOPQR'
valid_combinations = []

for n1, n2 in product(digits, digits):
    if n1 == n2: continue  # numbers can't repeat
    for l1, l2 in product(letters, letters):
        if l1 == l2: continue  # letters can't repeat
        guess = n1 + n2 + l1 + l2
        if check_conditions(guess):
            valid_combinations.append([n1, n2, l1, l2])

print(valid_combinations)