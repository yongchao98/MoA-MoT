import string

def check_password(nums, letters):
    # Convert numbers to strings for consistent handling
    n1, n2 = str(nums[0]), str(nums[1])
    l1, l2 = letters[0], letters[1]
    candidate = n1 + n2 + l1 + l2
    
    # Condition 1: 86WU
    # Both numbers wrong, one letter correct but wrong position, one letter wrong
    c1_nums_wrong = (n1 not in '86' and n2 not in '86')
    c1_letters = sum(l in 'WU' for l in letters) == 1
    if not (c1_nums_wrong and c1_letters):
        return False
    
    # Condition 2: 94MF
    # One number correct but wrong position, one wrong
    c2_nums = sum(n in '94' for n in [n1, n2]) == 1
    # Letters must be after M and F in alphabet
    c2_letters = all(l > 'M' for l in letters)
    if not (c2_nums and c2_letters):
        return False
    
    # Condition 3: 80UV
    # Both numbers and letters wrong
    c3_nums = (n1 not in '80' and n2 not in '80')
    c3_letters = ('U' not in letters and 'V' not in letters)
    if not (c3_nums and c3_letters):
        return False
    
    # Condition 4: 47TH
    # Both numbers wrong, one letter correct in correct position
    c4_nums = (n1 not in '47' and n2 not in '47')
    c4_letters = ('T' in letters and letters.index('T') == 0) or ('H' in letters and letters.index('H') == 1)
    if not (c4_nums and c4_letters):
        return False
    
    # Condition 5: 03SU
    # One number correct but wrong position, one wrong
    c5_nums = sum(n in '03' for n in [n1, n2]) == 1
    c5_numbers_size = int(n1) > 0 and int(n2) > 3
    c5_letters = ('S' not in letters and 'U' not in letters)
    if not (c5_nums and c5_numbers_size and c5_letters):
        return False
    
    # Additional check for number 9 from condition 2
    if '9' in [n1, n2]:
        if [n1, n2].index('9') == 0:  # 9 should be in second position
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