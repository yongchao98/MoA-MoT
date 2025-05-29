from itertools import product

def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    num = int(numbers)
    
    # Check number conditions
    if num <= 2: return False  # too small (from guess 5)
    if num >= 97: return False  # too large (from guess 1)
    
    # Known incorrect numbers
    wrong_numbers = ['97','70','17','72','02','54','40','09']
    if numbers in wrong_numbers: return False
    
    # From guess 6 (96CL): one number is correct and in position
    # Since 96 is too large, 6 must be in second position
    if numbers[1] != '6': return False
    
    # From guess 11 (09XR): one letter is correct and in position
    # R is too late in alphabet, so X must be the correct letter
    if letters[0] != 'X': return False
    
    # Check letter conditions
    wrong_letters = 'YSGFXPCLEW'  # can't use X twice
    if letters[1] in wrong_letters: return False
    
    # Letters can't be too early (before G) or too late (after R)
    if letters[1] < 'G' or letters[1] > 'R': return False
    
    # Additional number checks from feedback
    if num < 30: return False  # based on multiple "too small" feedbacks
    if num > 70: return False  # based on multiple "too large" feedbacks
    
    return True

# Generate all possible combinations
digits = '3456'  # refined based on constraints
letters = 'HIJKLMNOPQ'  # between G and R, excluding known wrong letters
valid_combinations = []

for n1 in digits:
    n2 = '6'  # we know second digit must be 6
    for l2 in letters:
        guess = [n1, n2, 'X', l2]
        if check_conditions(''.join(guess)):
            valid_combinations.append(guess)

print(valid_combinations)