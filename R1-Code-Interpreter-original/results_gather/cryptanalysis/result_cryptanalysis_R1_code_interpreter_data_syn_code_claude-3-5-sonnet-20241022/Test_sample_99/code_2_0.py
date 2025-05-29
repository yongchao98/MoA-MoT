from itertools import product

def check_feedback(guess, actual):
    n1, n2, l1, l2 = guess
    a1, a2, al1, al2 = actual
    
    # Check numbers
    correct_num_pos = sum(1 for i, j in zip([n1,n2], [a1,a2]) if i == j)
    correct_num = len(set([n1,n2]) & set([a1,a2]))
    
    # Check letters
    correct_let_pos = sum(1 for i, j in zip([l1,l2], [al1,al2]) if i == j)
    
    return correct_num_pos, correct_num, correct_let_pos

def is_too_late(letter):
    return letter >= 'T'

def is_too_early(letter):
    return letter <= 'G'

def validate_guess(candidate):
    n1, n2, l1, l2 = candidate
    
    # Condition 1: 87VX
    if n1 in '87' or n2 in '87':
        return False
    if not (is_too_late(l1) and is_too_late(l2)):  # This was wrong in original guess
        return True
    
    # Condition 2: 26AP
    if n1 in '26' or n2 in '26' or l1 in 'AP' or l2 in 'AP':
        return False
    
    # Condition 3: 98ZT - one number correct position
    nums_98 = check_feedback(['9','8','Z','T'], candidate)[0]
    if nums_98 != 1:
        return False
    
    # Condition 4: 98HL - one number correct position, one letter correct position
    nums_98hl, _, lets_98hl = check_feedback(['9','8','H','L'], candidate)
    if nums_98hl != 1 or lets_98hl != 1:
        return False
    
    # Condition 5: 93BN - one number correct position
    nums_93 = check_feedback(['9','3','B','N'], candidate)[0]
    if nums_93 != 1:
        return False
    
    # Condition 6: 37CG
    if n1 in '37' or n2 in '37':
        return False
    if not (is_too_early(l1) and is_too_early(l2)):  # This was wrong in original guess
        return True
    
    # Condition 7: 17KW
    if n1 in '17' or n2 in '17':
        return False
    
    # Condition 8: 92XC - one number correct position
    nums_92 = check_feedback(['9','2','X','C'], candidate)[0]
    if nums_92 != 1:
        return False
    
    # Condition 9: 30NX - both numbers too small
    if int(n1) <= 3 or int(n2) <= 0:
        return False
    
    # Condition 10: 74DJ - one number correct position
    nums_74 = check_feedback(['7','4','D','J'], candidate)[0]
    if nums_74 != 1:
        return False
    
    return True

# Generate all possible combinations
solutions = []
for n1, n2 in product('0123456789', repeat=2):
    if n1 == n2:
        continue
    for l1, l2 in product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2):
        if l1 == l2:
            continue
        candidate = [n1, n2, l1, l2]
        if validate_guess(candidate):
            solutions.append(candidate)

print(solutions)