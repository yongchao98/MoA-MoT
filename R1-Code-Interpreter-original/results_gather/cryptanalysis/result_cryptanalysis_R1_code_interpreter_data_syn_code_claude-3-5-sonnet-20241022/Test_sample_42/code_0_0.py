def check_conditions(nums, letters):
    # Convert to list for easier position checking
    guess = [str(nums[0]), str(nums[1]), letters[0], letters[1]]
    
    # Condition 4: 28JC - one number correct but wrong position, 8 too large
    if 8 in nums: return False
    if 2 not in nums: return False
    if guess.index(str(2)) == 0: return False  # 2 must be in second position
    
    # Condition 5: 98XJ - X in wrong position, both numbers too large
    if 9 in nums or 8 in nums: return False
    if 'X' not in letters: return False
    
    # Condition 6: 37XW - both X and W must be present but in opposite positions
    if 'X' not in letters or 'W' not in letters: return False
    
    # Letters that are too early in alphabet
    early_letters = set('ABCDEFGHIJKLMNOP')
    if any(l in early_letters for l in letters): return False
    
    # Numbers that are known to be incorrect
    wrong_nums = {0,3,4,5,7,9}
    if any(n in wrong_nums for n in nums): return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2: continue
        for l1 in 'QRSTUVWXYZ':
            for l2 in 'QRSTUVWXYZ':
                if l1 == l2: continue
                if check_conditions([n1,n2], [l1,l2]):
                    valid_solutions.append([n1,n2,l1,l2])

print(valid_solutions)