def check_conditions(guess):
    # Define all conditions
    conditions = [
        # 61OJ - both numbers wrong, letters too early
        lambda x: x[0] != '6' and x[1] != '1' and x[2] > 'O' and x[3] > 'J',
        
        # 98QN - one number correct in position, one too large, letters too early
        lambda x: ((x[0] == '9' and x[1] != '8') or (x[0] != '9' and x[1] == '8')) and x[2] > 'Q' and x[3] > 'N',
        
        # 58FC - one number correct in position, letters too early
        lambda x: ((x[0] == '5' and x[1] != '8') or (x[0] != '5' and x[1] == '8')),
        
        # 72YK - all incorrect
        lambda x: x[0] != '7' and x[1] != '2' and x[2] != 'Y' and x[3] != 'K',
        
        # 02ER - numbers too small, letters too early
        lambda x: int(x[0]) > 0 and int(x[1]) > 2 and x[2] > 'E' and x[3] > 'R',
        
        # 35QW - one number correct in position
        lambda x: ((x[0] == '3' and x[1] != '5') or (x[0] != '3' and x[1] == '5')),
        
        # 71KJ - all incorrect, letters too early
        lambda x: x[0] != '7' and x[1] != '1' and x[2] > 'K' and x[3] > 'J',
        
        # 79PX - numbers wrong, one letter correct in position
        lambda x: x[0] != '7' and x[1] != '9' and (x[2] == 'P' or x[3] == 'X'),
        
        # 86TV - one number correct wrong position, one letter correct in position
        lambda x: (x[1] == '8' or x[0] == '8') and (x[2] == 'T' or x[3] == 'V')
    ]
    
    # Check if all conditions are met
    return all(cond(guess) for cond in conditions)

# Generate possible combinations
def find_password():
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:  # Skip if numbers repeat
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:  # Skip if letters repeat
                        continue
                    guess = n1 + n2 + l1 + l2
                    if check_conditions(guess):
                        print(f"Found solution: {[n1, n2, l1, l2]}")

find_password()