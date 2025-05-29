from itertools import permutations

def check_guess(guess):
    n1, n2, l1, l2 = guess
    
    # Convert numbers to integers for comparison
    num1, num2 = int(n1), int(n2)
    
    # Condition 1: 92BJ - both numbers and letters wrong
    if n1 in '92' or n2 in '92' or l1 in 'BJ' or l2 in 'BJ':
        return False
    
    # Condition 2: 06RU - one number correct in position
    c2_num_match = (n1 == '0' or n2 == '6')
    if not c2_num_match:
        return False
    if l1 >= 'R' or l2 >= 'R':  # letters too late
        return False
        
    # Condition 3: 84CV - both numbers wrong, one letter correct in position
    if n1 in '84' or n2 in '84':
        return False
    if not (l1 == 'C' or l2 == 'C'):  # one letter must be C
        return False
        
    # Condition 4: 75VY - one number correct in position
    c4_num_match = (n1 == '7' or n2 == '5')
    if not c4_num_match:
        return False
    if l1 in 'VY' or l2 in 'VY':
        return False
        
    # Condition 5: 16UD - one number correct in position
    if l1 in 'UD' or l2 in 'UD':
        return False
        
    # Condition 6: 98NG - numbers too large
    if num1 >= 9 or num2 >= 8:
        return False
    if l1 in 'NG' or l2 in 'NG':
        return False
        
    # Condition 7: 03QD - numbers too small, Q must be present
    if num1 <= 0 or num2 <= 3:
        return False
    if 'Q' not in [l1, l2]:  # Q must be present
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

valid_solutions = []
for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        guess = [n1, n2, l1, l2]
        if check_guess(guess):
            valid_solutions.append(guess)

print(valid_solutions)