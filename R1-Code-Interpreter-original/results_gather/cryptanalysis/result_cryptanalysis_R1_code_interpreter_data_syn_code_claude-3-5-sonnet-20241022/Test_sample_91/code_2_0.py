def check_solution(guess):
    n1, n2, l1, l2 = guess
    
    # Check if numbers or letters repeat
    if n1 == n2 or l1 == l2:
        return False
    
    # Check numbers are less than 4
    if n1 >= 4 or n2 >= 4:
        return False
    
    # Check letters are only I and M
    if set([l1, l2]) != set(['I', 'M']):
        return False
    
    # Condition 2 (10FX): if 0 is present, it must be in wrong position compared to guess
    if 0 in [n1, n2]:
        if n2 == 0:  # 0 is in second position in guess 2
            return False
    zero_check = sum(n in [1,0] for n in [n1,n2]) == 1
    if not zero_check:
        return False
    
    # Condition 3 (50ZK): one of 0,5 must be in wrong position
    zero_five_check = sum(n in [0] for n in [n1,n2]) == 1  # we already know numbers < 4
    if not zero_five_check:
        return False
    
    # Condition 6 (39KL): one of 3,9 must be in wrong position
    three_check = sum(n in [3] for n in [n1,n2]) == 1  # we already know numbers < 4
    if not three_check:
        return False
    
    # From guess 4 (48IM): I and M must be in wrong positions
    if l1 == 'I' and l2 == 'M':  # Same positions as guess 4
        return False
    
    # From guess 5 (86NM): M cannot be in last position
    if l2 == 'M':  # M in same position as guess 5
        return False
    
    return True

# Try all possible combinations
valid_solutions = []
for n1 in range(4):
    for n2 in range(4):
        for l1 in ['I', 'M']:
            for l2 in ['I', 'M']:
                guess = [n1, n2, l1, l2]
                if check_solution(guess):
                    valid_solutions.append(guess)
print(valid_solutions)