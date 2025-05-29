def check_conditions(numbers, letters):
    guess = numbers + letters
    
    # Helper function to count correct positions and values
    def check_numbers(test_nums, actual_nums):
        correct_pos = sum(1 for i in range(2) if test_nums[i] == actual_nums[i])
        correct_val = sum(1 for x in test_nums if x in actual_nums)
        return correct_pos, correct_val
    
    # Condition 1: 52OL - both numbers incorrect, letters too early
    if '5' in numbers or '2' in numbers:
        return False
    
    # Condition 2: 20FP - one number correct position
    n_pos, n_val = check_numbers(['2', '0'], numbers)
    if n_pos != 1:
        return False
    
    # Condition 3: 08EX - one number correct but wrong position
    n_pos, n_val = check_numbers(['0', '8'], numbers)
    if n_val != 1 or n_pos != 0:
        return False
    
    # Condition 4: 58RN - both numbers too large
    if not (int(numbers[0]) < 5 and int(numbers[1]) < 5):
        return False
    
    # Condition 5: 64UM - one number wrong pos, one letter correct pos
    if 'U' not in letters:
        return False
    n_pos, n_val = check_numbers(['6', '4'], numbers)
    if n_val != 1 or n_pos != 0:
        return False
    
    # Condition 6: 02PV - one number wrong pos
    n_pos, n_val = check_numbers(['0', '2'], numbers)
    if n_val != 1 or n_pos != 0:
        return False
    
    # Condition 7: 62CO - both numbers incorrect
    if '6' in numbers or '2' in numbers:
        return False
    
    # Condition 8: 91CD - both numbers incorrect
    if '9' in numbers or '1' in numbers:
        return False
    
    # Condition 9: 96SR - numbers too large, S must be present but wrong pos
    if 'S' not in letters:
        return False
    
    return True

# Generate and test all possible combinations
solutions = []
for d1 in range(5):  # 0-4 based on "too large" conditions
    for d2 in range(5):
        if d1 == d2:
            continue
        for l1 in ['U', 'S']:  # We know U and S must be present
            for l2 in ['U', 'S']:
                if l1 == l2:
                    continue
                numbers = [str(d1), str(d2)]
                letters = [l1, l2]
                if check_conditions(numbers, letters):
                    solutions.append([str(d1), str(d2), l1, l2])

print(solutions)