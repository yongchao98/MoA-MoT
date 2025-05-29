from itertools import permutations

def check_conditions(nums, letters):
    # Convert to strings for easier handling
    n1, n2 = str(nums[0]), str(nums[1])
    l1, l2 = letters[0], letters[1]
    
    # Condition checks
    if any(x in [6,2,7,1,3,8,9] for x in nums):
        return False
    
    # Numbers can't be both too small (from clue 3)
    if all(x < 4 for x in nums):
        return False
        
    # Numbers can't be both too large (from clue 4)
    if all(x > 5 for x in nums):
        return False
    
    # From clue 7: one number must be 3 or 5, but not both
    if not ((nums[0] == 3 and nums[1] != 5) or (nums[0] == 5 and nums[1] != 3) or
            (nums[1] == 3 and nums[0] != 5) or (nums[1] == 5 and nums[0] != 3)):
        return False
    
    # Letter conditions
    if any(x in 'SNBKWRTMDU' for x in letters):
        return False
    
    # A and G must be present but not in same positions as in '67AG'
    if 'A' not in letters or 'G' not in letters:
        return False
    if letters == ('A', 'G'):
        return False
    
    return True

# Generate possible numbers (0-9)
possible_nums = list(range(10))
num_pairs = list(permutations(possible_nums, 2))

# Generate possible letters (A-Z)
possible_letters = [chr(i) for i in range(65, 91)]  # A-Z
letter_pairs = list(permutations(['A', 'G'], 2))

# Test all combinations
solutions = []
for nums in num_pairs:
    for letters in letter_pairs:
        if check_conditions(nums, letters):
            solutions.append([str(nums[0]), str(nums[1]), letters[0], letters[1]])

print(solutions)