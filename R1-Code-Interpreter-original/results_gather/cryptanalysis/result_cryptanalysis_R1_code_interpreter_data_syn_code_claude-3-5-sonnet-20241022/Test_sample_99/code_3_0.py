def validate_combination(n1, n2, l1, l2):
    # Helper function to check if a letter is too late in alphabet
    def is_too_late(letter):
        return letter in 'TUVWXYZ'
    
    # Helper function to check if a letter is too early in alphabet
    def is_too_early(letter):
        return letter in 'ABCDEFG'
    
    # 1. 87VX: both numbers wrong, both letters too late
    if n1 in '87' or n2 in '87':
        return False
    
    # 2. 26AP: both numbers and letters wrong
    if n1 in '26' or n2 in '26' or l1 in 'AP' or l2 in 'AP':
        return False
    
    # 3. 98ZT: one number correct and in position
    correct_pos_98 = (n1 == '9') + (n2 == '8')
    if correct_pos_98 != 1:
        return False
    
    # 4. 98HL: one number correct position, one letter correct position
    correct_pos_98hl = (n1 == '9') + (n2 == '8')
    correct_pos_hl = (l1 == 'H') + (l2 == 'L')
    if correct_pos_98hl != 1 or correct_pos_hl != 1:
        return False
    
    # 5. 93BN: one number correct position, one too small
    if not ((n1 == '9' and int(n2) > 3) or (n2 == '3' and int(n1) > 3)):
        return False
    
    # 6. 37CG: both numbers wrong, letters too early
    if n1 in '37' or n2 in '37':
        return False
    
    # 7. 17KW: both numbers wrong
    if n1 in '17' or n2 in '17':
        return False
    
    # 8. 92XC: one number correct position, one too small
    correct_pos_92 = (n1 == '9') + (n2 == '2')
    if correct_pos_92 != 1:
        return False
    
    # 9. 30NX: both numbers too small
    if int(n1) <= 3 or int(n2) <= 0:
        return False
    
    # 10. 74DJ: one number correct position
    correct_pos_74 = (n1 == '7') + (n2 == '4')
    if correct_pos_74 != 1:
        return False
    
    # Additional letter constraints from all conditions
    if is_too_late(l1) or is_too_late(l2):
        return False
    if is_too_early(l1) or is_too_early(l2):
        return False
    
    return True

# Test all possible combinations
solutions = []
for n1 in '0123456789':
    for n2 in '0123456789':
        if n1 == n2:
            continue
        for l1 in 'HIJKLMNOPQRS':  # Middle alphabet letters
            for l2 in 'HIJKLMNOPQRS':
                if l1 == l2:
                    continue
                if validate_combination(n1, n2, l1, l2):
                    solutions.append([n1, n2, l1, l2])

print(solutions)