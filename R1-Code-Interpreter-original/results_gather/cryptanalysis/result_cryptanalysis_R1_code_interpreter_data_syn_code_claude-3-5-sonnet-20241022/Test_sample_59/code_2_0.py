def check_all_conditions(n1, n2, l1, l2):
    # Basic conditions
    if l2 != 'O': return False
    if n1 != '1': return False
    if l1 < 'H' or l1 > 'T': return False
    if l1 == 'O': return False  # Letters can't repeat
    
    # From 09YT: one number correct but wrong position
    # This means either 0 or 9 should be in the correct position
    # Since neither is correct, these numbers are eliminated
    if n2 in '09': return False
    
    # From 15LD: one number correct but wrong position
    # We know 1 is correct in first position, so 5 can't be correct
    if n2 == '5': return False
    
    # From 78LK, 48TP: both numbers incorrect
    if n2 in '478': return False
    
    # From 65RX: both numbers incorrect
    if n2 in '65': return False
    
    # From 03CO and 72NO: O is correct in position and other letter is too early
    # This means l1 must be after C and N in alphabet
    if l1 <= 'N': return False
    
    # From 21XS: letters are too late in alphabet
    if l1 >= 'S': return False
    
    # From all guesses where letters are "too early"
    if l1 <= 'L': return False  # L is mentioned multiple times as too early
    
    # From all guesses where letters are "too late"
    if l1 >= 'T': return False  # T is mentioned as too late
    
    return True

valid_combinations = []
first_num = '1'
possible_second_nums = '23'  # After eliminating all impossible numbers
possible_first_letters = 'MNOPQR'  # After considering alphabet position constraints

for n2 in possible_second_nums:
    for l1 in possible_first_letters:
        if check_all_conditions(first_num, n2, l1, 'O'):
            valid_combinations.append([first_num, n2, l1, 'O'])

print(valid_combinations)