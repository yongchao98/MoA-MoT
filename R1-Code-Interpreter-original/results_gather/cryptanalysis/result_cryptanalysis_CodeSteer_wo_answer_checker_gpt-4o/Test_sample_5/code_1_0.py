from itertools import permutations

def is_valid_combination(numbers, letters, combination):
    num1, num2, let1, let2 = combination
    
    # Condition 1: 20KJ
    if num1 <= 2 or num2 <= 0 or let1 in "KJ" or let2 in "KJ":
        return False
    
    # Condition 2: 60SC
    if not ((num1 == 6 and num2 != 0) or (num2 == 6 and num1 != 0)):
        return False
    if not ((let1 == 'S' and let2 != 'C') or (let2 == 'S' and let1 != 'C')):
        return False
    
    # Condition 3: 39OQ
    if not ((num1 == 3 and num2 < 9) or (num2 == 3 and num1 < 9)):
        return False
    if let1 in "OQ" or let2 in "OQ":
        return False
    
    # Condition 4: 09BV
    if num1 in [0, 9] or num2 in [0, 9]:
        return False
    if not ((let1 == 'B' and let2 != 'V') or (let2 == 'B' and let1 != 'V')):
        return False
    
    return True

def find_password():
    numbers = range(10)
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    for num_comb in permutations(numbers, 2):
        for let_comb in permutations(letters, 2):
            combination = (*num_comb, *let_comb)
            if is_valid_combination(numbers, letters, combination):
                return combination

password = find_password()
print(f"<<< {list(map(str, password))} >>>")