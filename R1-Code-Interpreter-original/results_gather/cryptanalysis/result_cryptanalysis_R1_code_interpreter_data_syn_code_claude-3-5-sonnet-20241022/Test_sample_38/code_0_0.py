def check_conditions(nums, letters):
    # Convert guess to list for easier comparison
    guess = [str(nums[0]), str(nums[1]), letters[0], letters[1]]
    
    # Check condition 1: 47YQ
    if not (int(nums[0]) < 4 and int(nums[1]) < 4):
        return False
    
    # Check condition 2: 91LQ
    # One letter must be correct and in position
    c2_letters = [letters[i] == 'L' and i == 0 or letters[i] == 'Q' and i == 1 for i in range(2)]
    if sum(c2_letters) != 1:
        return False
    if nums[0] == '9' or nums[1] == '1':
        return False
    
    # Check condition 3: 23TP
    # One number must be correct but in wrong position
    c3_nums = ['2', '3']
    if sum(n in c3_nums for n in nums) != 1:
        return False
    if any(nums.index(n) == c3_nums.index(n) for n in c3_nums if n in nums):
        return False
    
    # Check condition 4: 09KG
    # One number correct but wrong position, one number too large
    c4_nums = ['0', '9']
    if sum(n in c4_nums for n in nums) != 1:
        return False
    if '9' in nums:  # 9 is too large
        return False
    if any(nums.index(n) == c4_nums.index(n) for n in c4_nums if n in nums):
        return False
    
    # One letter correct and in position
    c4_letters = [letters[i] == 'K' and i == 0 or letters[i] == 'G' and i == 1 for i in range(2)]
    if sum(c4_letters) != 1:
        return False
    
    return True

# Generate all possible combinations
possible_nums = []
for i in range(4):
    for j in range(4):
        if i != j:
            possible_nums.append((i, j))

possible_letters = []
for l1 in 'ABCDEFGHIJKL':  # Letters before M (as later ones are too late)
    for l2 in 'ABCDEFGHIJKL':
        if l1 != l2:
            possible_letters.append((l1, l2))

# Find the solution
for nums in possible_nums:
    for letters in possible_letters:
        if check_conditions(nums, letters):
            print([str(nums[0]), str(nums[1]), letters[0], letters[1]])