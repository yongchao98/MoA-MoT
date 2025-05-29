def verify_final_order(d1, d2, l1, l2):
    # Condition 8: 58IJ - one letter correct but in the wrong position
    # This is crucial - if I is in second position (l2), 
    # then it would be in the same position as in 58IJ, which contradicts "wrong position"
    # Therefore, I must be in first letter position (l1)
    
    # Condition 3: 92FQ - both letters incorrect
    # If F is in the same position as in our answer, it would be "correct but wrong position"
    # Since it's "incorrect", F must be in different position than in 92FQ
    # In 92FQ, F is in first position, so in our answer F must be second position
    
    if l1 != 'I' or l2 != 'F':
        return False
    
    return True

valid_combinations = []
d1, d2 = 3, 0
for l1 in ['F', 'I']:
    for l2 in ['F', 'I']:
        if l1 == l2:
            continue
        if verify_final_order(d1, d2, l1, l2):
            valid_combinations.append([d1, d2, l1, l2])

print("Final combination:")
for combo in valid_combinations:
    print(combo)