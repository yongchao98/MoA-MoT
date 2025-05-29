def verify_all_conditions(d1, d2, l1, l2):
    # Condition 1: 16QS - both numbers incorrect; both letters incorrect and too late
    if d1 in [1, 6] or d2 in [1, 6]:
        return False
    
    # Condition 2: 23YQ - one number correct but wrong position
    # Since we know d1=3, this means 2 cannot be d2
    if d2 == 2:
        return False
    
    # Condition 3: 92FQ - both numbers incorrect
    # This confirms F could be one of our letters
    
    # Condition 4: 13MS - one number correct but wrong position
    # Since d1=3, this means 1 cannot be d2
    if d2 == 1:
        return False
    
    # Condition 8: 58IJ - one letter correct but wrong position
    # We already know one letter must be I
    
    # Condition 13: 90QC - one number correct and in correct position
    # Since both our numbers are 3,0, and 9 is too large, 0 must be correct
    if d2 != 0:
        return False
    
    # From multiple conditions, we know:
    # - Letters can't be too late in alphabet (after P)
    # - One letter must be I
    # - The other letter must be F (from condition 3 and other constraints)
    if not ((l1 == 'I' and l2 == 'F') or (l1 == 'F' and l2 == 'I')):
        return False
    
    return True

valid_combinations = []
# We know first digit is 3 and second is 0
d1, d2 = 3, 0
for l1 in [chr(i) for i in range(ord('D'), ord('P')+1)]:
    for l2 in [chr(i) for i in range(ord('D'), ord('P')+1)]:
        if l1 == l2:
            continue
        if verify_all_conditions(d1, d2, l1, l2):
            valid_combinations.append([d1, d2, l1, l2])

print("Final possible combinations:")
for combo in valid_combinations:
    print(combo)