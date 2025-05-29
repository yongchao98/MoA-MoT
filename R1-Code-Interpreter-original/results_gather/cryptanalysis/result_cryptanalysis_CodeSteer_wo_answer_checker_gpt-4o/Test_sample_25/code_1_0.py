from itertools import permutations

# Possible numbers and letters
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Feedback conditions
feedback = [
    ('91SF', (0, 1)),  # Both numbers incorrect, one letter correct and in position
    ('03HF', (1, 0)),  # One number correct and in position, one number too small
    ('75CP', (0, 1)),  # Both numbers incorrect, one letter correct but wrong position
    ('64YB', (1, 0))   # One number correct and in position, one number incorrect
]

def is_valid_combination(num1, num2, let1, let2):
    # Check against each feedback
    for guess, (num_feedback, let_feedback) in feedback:
        guess_nums = guess[:2]
        guess_lets = guess[2:]
        
        # Check numbers
        num_correct = sum((num1 == guess_nums[0], num2 == guess_nums[1]))
        if num_correct != num_feedback:
            return False
        
        # Check letters
        let_correct = sum((let1 == guess_lets[0], let2 == guess_lets[1]))
        if let_correct != let_feedback:
            return False
        
        # Additional checks for too small or wrong position
        if guess == '03HF' and num1 != '0':
            return False
        if guess == '75CP' and (let1 == 'C' or let2 == 'P'):
            return False
        if guess == '64YB' and num2 == '4':
            return False
    
    return True

# Iterate through all possible combinations
for num1, num2 in permutations(numbers, 2):
    for let1, let2 in permutations(letters, 2):
        if is_valid_combination(num1, num2, let1, let2):
            password = [num1, num2, let1, let2]
            print(f"<<< {password} >>>")
            break