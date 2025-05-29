import string

def check_password(nums, letters):
    # Convert numbers to strings for consistent handling
    n1, n2 = str(nums[0]), str(nums[1])
    l1, l2 = letters[0], letters[1]
    
    # Condition 1: 86WU
    # Both numbers wrong, one letter correct but wrong position, one letter wrong
    if n1 in '86' or n2 in '86':
        return False
    if sum(l in 'WU' for l in letters) != 1:
        return False
    if 'W' in letters and letters.index('W') == 1:  # W must be in wrong position
        return False
    if 'U' in letters and letters.index('U') == 0:  # U must be in wrong position
        return False
    
    # Condition 2: 94MF
    # One number correct but wrong position, one wrong, both letters too early
    if sum(n in '94' for n in [n1, n2]) != 1:
        return False
    if '9' in [n1, n2] and [n1, n2].index('9') == 0:  # If 9 present, must be in wrong position
        return False
    if '4' in [n1, n2] and [n1, n2].index('4') == 1:  # If 4 present, must be in wrong position
        return False
    if any(l <= 'M' for l in letters):  # Letters must be after M
        return False
    
    # Condition 3: 80UV
    # Both numbers and letters wrong
    if n1 in '80' or n2 in '80':
        return False
    if 'U' in letters or 'V' in letters:
        return False
    
    # Condition 4: 47TH
    # Both numbers wrong, one letter correct and in position, one letter too early
    if n1 in '47' or n2 in '47':
        return False
    if sum(l in 'TH' for l in letters) != 1:
        return False
    if 'T' in letters and letters.index('T') != 0:  # T must be in first position if present
        return False
    if 'H' in letters:  # H is too early in alphabet
        return False
    
    # Condition 5: 03SU
    # One number correct but wrong position, one too small, both letters wrong
    if sum(n in '03' for n in [n1, n2]) != 1:
        return False
    if int(n1) <= 0 or int(n2) <= 3:  # Numbers must be greater than these values
        return False
    if 'S' in letters or 'U' in letters:
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:  # No repeating numbers
            for l1 in string.ascii_uppercase:
                for l2 in string.ascii_uppercase:
                    if l1 != l2:  # No repeating letters
                        if check_password([n1, n2], [l1, l2]):
                            valid_solutions.append([str(n1), str(n2), l1, l2])

print(valid_solutions)