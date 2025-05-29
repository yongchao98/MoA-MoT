# Deduced numbers and letters
numbers = [2, 4]
letters = ['W', 'R']

# Function to check the conditions
def check_conditions():
    conditions = [
        (86, 'AF', 'both numbers are incorrect; both letters are incorrect'),
        (98, 'LF', 'both numbers are incorrect and too large; both letters are incorrect'),
        (20, 'XK', 'one number is correct and in the correct position; one number is incorrect and too small; both letters are incorrect'),
        (38, 'TM', 'both numbers are incorrect; both letters are incorrect'),
        (81, 'JK', 'both numbers are incorrect; both letters are incorrect'),
        (23, 'VQ', 'one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect'),
        (85, 'IC', 'both numbers are incorrect; both letters are incorrect'),
        (03, 'OE', 'both numbers are incorrect; both letters are incorrect'),
        (78, 'WK', 'one number is correct but in the wrong position; one number is incorrect and too large; one letter is correct and in the correct position; one letter is incorrect and incorrect'),
        (68, 'KC', 'both numbers are incorrect; both letters are incorrect'),
        (47, 'OM', 'one number is correct and in the correct position; one number is incorrect and incorrect; both letters are incorrect'),
        (58, 'ZQ', 'both numbers are incorrect; both letters are incorrect'),
        (49, 'DW', 'both numbers are incorrect; one letter is correct but in the wrong position; one letter is incorrect and incorrect')
    ]
    
    for num, let, feedback in conditions:
        num_str = str(num)
        let_str = let
        num_correct = sum(1 for i, n in enumerate(num_str) if int(n) in numbers and numbers[i] == int(n))
        num_wrong_pos = sum(1 for n in num_str if int(n) in numbers) - num_correct
        num_incorrect = 2 - num_correct - num_wrong_pos
        
        let_correct = sum(1 for i, l in enumerate(let_str) if l in letters and letters[i] == l)
        let_wrong_pos = sum(1 for l in let_str if l in letters) - let_correct
        let_incorrect = 2 - let_correct - let_wrong_pos
        
        if 'both numbers are incorrect' in feedback and num_incorrect != 2:
            return False
        if 'both letters are incorrect' in feedback and let_incorrect != 2:
            return False
        if 'one number is correct and in the correct position' in feedback and num_correct != 1:
            return False
        if 'one number is incorrect and too small' in feedback and num_wrong_pos != 1:
            return False
        if 'one number is correct but in the wrong position' in feedback and num_wrong_pos != 1:
            return False
        if 'one letter is correct and in the correct position' in feedback and let_correct != 1:
            return False
        if 'one letter is correct but in the wrong position' in feedback and let_wrong_pos != 1:
            return False
    
    return True

# Check if the deduced numbers and letters satisfy all conditions
print(check_conditions())