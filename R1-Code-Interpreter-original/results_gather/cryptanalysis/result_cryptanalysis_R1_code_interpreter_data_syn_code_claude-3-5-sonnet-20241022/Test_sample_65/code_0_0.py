from itertools import permutations

def check_conditions(password):
    n1, n2, l1, l2 = password
    
    # Convert numbers to strings for comparison
    n1, n2 = str(n1), str(n2)
    
    # Condition 1 & 3: 95 are correct numbers but wrong positions
    if not (n1 in '95' and n2 in '95'):
        return False
    
    # Condition 7: one number correct in position
    # Since we know numbers are 9 and 5, and one is in correct position
    # Testing guess "19SZ" - one number correct in position
    if (n1 == '1' or n2 == '9'):
        return False
    
    # Condition 6: 16OY - O is correct in position
    if l1 != 'O':
        return False
    
    # Condition 41FV - V is correct but wrong position
    if l2 != 'V':
        return False
    
    return True

# Generate all possible combinations
numbers = [5, 9]  # We know these are the numbers
letters = ['O', 'V']  # We know O and V are the letters

all_combinations = list(permutations(numbers, 2))

for nums in all_combinations:
    password = [nums[0], nums[1], 'O', 'V']
    if check_conditions(password):
        print(password)