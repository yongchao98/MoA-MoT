def verify_solution(nums, lets):
    # Guess 1: 80KR - one number correct but wrong pos, one too small, both letters wrong
    if nums.count('8') + nums.count('0') != 1:  # exactly one of 8,0 should be correct
        return False
    if '8' in nums and nums[0] == '8':  # 8 should be in wrong position
        return False
    if 'K' in lets or 'R' in lets:  # K and R should not be present
        return False
    
    # Guess 2: 87NP - one number correct but wrong pos, one incorrect, one letter correct in pos
    if nums.count('8') + nums.count('7') != 1:  # exactly one of 8,7 should be correct
        return False
    if '8' in nums and nums[0] == '8':  # 8 should be in wrong position
        return False
    if (lets[0] != 'N' and lets[1] != 'P'):  # either N or P must be in correct position
        return False
    
    # Guess 3: 47AZ - both numbers wrong, one letter correct in pos, Z too late
    if '4' in nums or '7' in nums:
        return False
    if lets[0] != 'A':  # A must be in first position (since Z is too late)
        return False
    
    # Guess 4: 15RF - both numbers too small, both letters wrong
    if int(nums[0]) <= 1 or int(nums[1]) <= 5:
        return False
    if 'R' in lets or 'F' in lets:
        return False
    
    return True

# Test all possible combinations
valid_solutions = []
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        nums = [str(n1), str(n2)]
        # We know A must be first letter and P must be second letter
        lets = ['A', 'P']
        if verify_solution(nums, lets):
            valid_solutions.append(nums + lets)

print(valid_solutions)