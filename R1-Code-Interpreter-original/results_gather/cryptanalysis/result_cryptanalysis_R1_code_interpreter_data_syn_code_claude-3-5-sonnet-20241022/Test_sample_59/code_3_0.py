def check_final_conditions(n1, n2, l1, l2):
    # Basic conditions from previous analysis
    if n1 != '1' or l2 != 'O': return False
    if l1 not in 'PQR': return False
    if n2 not in '23': return False
    
    # From 21XS: one number is correct and in correct position (1), 
    # one number is incorrect - this means 2 cannot be the second number
    if n2 == '2': return False
    
    # From 03CO and 72NO: 
    # O is correct in position and other letter is too early in alphabet
    # This suggests P is the optimal choice as it's right after N
    if l1 != 'P': return False
    
    # Additional verification:
    # - Q and R would be too late in alphabet (closer to S, T, X which are marked as too late)
    # - P is the perfect position being right after N (which is too early) 
    #   but well before S, T, X (which are too late)
    
    return True

valid_combinations = []
for n2 in '23':
    for l1 in 'PQR':
        if check_final_conditions('1', n2, l1, 'O'):
            valid_combinations.append(['1', n2, l1, 'O'])

print(valid_combinations)