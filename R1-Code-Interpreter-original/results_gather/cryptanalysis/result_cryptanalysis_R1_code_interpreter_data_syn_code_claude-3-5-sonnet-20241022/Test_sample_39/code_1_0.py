def validate_password(password):
    # Check first guess (87ID)
    n1, n2, l1, l2 = password
    
    # For 87ID
    # Check numbers
    correct_num_wrong_pos_87 = (n1 == '7' or n2 == '7') and '7' in "87"
    too_large_count_87 = sum(1 for x in "87" if int(x) > int(n1)) + sum(1 for x in "87" if int(x) > int(n2))
    
    # Check letters
    correct_letter_wrong_pos_87 = ((l1 in "ID" or l2 in "ID") and (l1 != 'I' or l2 != 'D'))
    too_early_letter_87 = sum(1 for x in "ID" if x < l1) + sum(1 for x in "ID" if x < l2) == 1
    
    # For 50FH
    # Check numbers
    too_small_50 = int(n1) > 5 and int(n2) > 0
    
    # Check letters
    correct_letter_wrong_pos_50 = ((l1 in "FH" or l2 in "FH") and (l1 != 'F' or l2 != 'H'))
    too_early_letter_50 = sum(1 for x in "FH" if x < l1) + sum(1 for x in "FH" if x < l2) == 1
    
    return (correct_num_wrong_pos_87 and too_large_count_87 == 1 and 
            correct_letter_wrong_pos_87 and too_early_letter_87 and
            too_small_50 and correct_letter_wrong_pos_50 and too_early_letter_50)

# Generate all possible combinations
valid_passwords = []
for n1 in range(1, 10):
    for n2 in range(1, 10):
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2:
                        password = [str(n1), str(n2), l1, l2]
                        if validate_password(password):
                            valid_passwords.append(password)

print(valid_passwords)