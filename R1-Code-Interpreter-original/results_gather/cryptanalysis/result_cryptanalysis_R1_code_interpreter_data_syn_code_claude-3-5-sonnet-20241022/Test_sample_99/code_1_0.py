from itertools import product

def is_valid_guess(candidate):
    num1, num2, let1, let2 = candidate
    
    # Helper function to check if a letter is "too late" in alphabet
    def is_too_late(letter):
        return letter in 'TUVWXYZ'
    
    # Helper function to check if a letter is "too early" in alphabet
    def is_too_early(letter):
        return letter in 'ABCDEFG'
    
    # Check all conditions
    
    # Condition 1: 87VX - both numbers wrong, both letters too late
    if num1 in '87' or num2 in '87' or not (is_too_late(let1) and is_too_late(let2)):
        return False
    
    # Condition 2: 26AP - both numbers and letters wrong
    if num1 in '26' or num2 in '26' or let1 in 'AP' or let2 in 'AP':
        return False
    
    # Condition 3: 98ZT - one number correct position, both letters too late
    correct_98 = (num1 == '9' and num2 != '8') or (num2 == '8' and num1 != '9')
    if not correct_98 or not (is_too_late(let1) and is_too_late(let2)):
        return False
    
    # Condition 4: 98HL - one number correct position, one letter correct position
    correct_98_2 = (num1 == '9' and num2 != '8') or (num2 == '8' and num1 != '9')
    correct_HL = (let1 == 'H' and let2 != 'L') or (let2 == 'L' and let1 != 'H')
    if not (correct_98_2 and correct_HL):
        return False
    
    # Condition 5: 93BN - one number correct position, other number too small
    if not ((num1 == '9' and int(num2) > 3) or (num2 == '3' and num1 != '9')):
        return False
    
    # Condition 6: 37CG - both numbers wrong, letters too early
    if num1 in '37' or num2 in '37' or not (is_too_early(let1) and is_too_early(let2)):
        return False
    
    # Condition 7: 17KW - both numbers wrong, letters too late
    if num1 in '17' or num2 in '17' or not (is_too_late(let1) and is_too_late(let2)):
        return False
    
    # Condition 8: 92XC - one number correct position, other too small
    if not ((num1 == '9' and int(num2) > 2) or (num2 == '2' and num1 != '9')):
        return False
    
    # Condition 9: 30NX - both numbers too small, letters too late
    if not (int(num1) > 3 and int(num2) > 0) or not (is_too_late(let1) and is_too_late(let2)):
        return False
    
    # Condition 10: 74DJ - one number correct position
    correct_74 = (num1 == '7' and num2 != '4') or (num2 == '4' and num1 != '7')
    if not correct_74:
        return False
    
    return True

# Generate and test all possible combinations
valid_solutions = []
for nums in product('0123456789', repeat=2):
    if nums[0] == nums[1]:  # Skip repeating numbers
        continue
    for lets in product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=2):
        if lets[0] == lets[1]:  # Skip repeating letters
            continue
        candidate = [nums[0], nums[1], lets[0], lets[1]]
        if is_valid_guess(candidate):
            valid_solutions.append(candidate)

print(valid_solutions)