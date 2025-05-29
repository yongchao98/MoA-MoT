def check_all_conditions(nums, letters):
    # Convert to list for position checking
    combo = [nums[0], nums[1], letters[0], letters[1]]
    
    # Basic constraints
    if nums[0] == nums[1] or letters[0] == letters[1]: return False
    
    # Numbers constraints
    wrong_nums = {0,3,4,5,7,8,9}
    if any(n in wrong_nums for n in nums): return False
    if 2 not in nums: return False  # 2 must be present
    
    # Letters constraints
    if set(letters) != {'X', 'W'}: return False  # Must be X and W
    
    # Position constraints from guess 4 (28JC)
    if combo.index(2) == 0: return False  # 2 can't be in first position
    
    # Position constraints from guess 5 (98XJ)
    x_pos_5 = 2  # X position in guess 5
    if letters.index('X') == x_pos_5: return False  # X can't be in same position as in guess 5
    
    # Position constraints from guess 6 (37XW)
    x_pos_6 = 2  # X position in guess 6
    w_pos_6 = 3  # W position in guess 6
    if letters.index('X') == x_pos_6 or letters.index('W') == w_pos_6: return False
    
    return True

# Generate all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2: continue
        # Try both orders of X and W
        for letters in [['X', 'W'], ['W', 'X']]:
            if check_all_conditions([n1, n2], letters):
                valid_solutions.append([n1, n2, letters[0], letters[1]])

print(valid_solutions)