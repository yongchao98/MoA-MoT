def check_all_conditions(nums, letters):
    guess = [str(nums[0]), str(nums[1]), letters[0], letters[1]]
    
    # Condition 4: 28JC - 2 must be present but in wrong position
    if 2 not in nums: return False
    if nums[1] == 2: return False  # 2 can't be in second position (from guess 4)
    
    # Condition 5: 98XJ - X in wrong position
    if 'X' not in letters: return False
    if letters[0] == 'X': return False  # X can't be in first letter position
    
    # Condition 6: 37XW - both X and W in wrong positions
    if 'W' not in letters: return False
    if letters[0] == 'W' or letters[1] == 'X': return False  # From positions in guess 6
    
    # Numbers that are known to be incorrect
    wrong_nums = {0,3,4,5,7,8,9}
    if any(n in wrong_nums for n in nums): return False
    
    # Letters that are too early in alphabet
    early_letters = set('ABCDEFGHIJKLMNOP')
    if any(l in early_letters for l in letters): return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2: continue
        for l1 in 'QRSTUVWXYZ':
            for l2 in 'QRSTUVWXYZ':
                if l1 == l2: continue
                if check_all_conditions([n1,n2], [l1,l2]):
                    valid_solutions.append([n1,n2,l1,l2])

print(valid_solutions)