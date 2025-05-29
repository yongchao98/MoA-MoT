# Possible numbers and letters based on the feedback
numbers = ['9', '5', '1', '2']
letters = ['W', 'C', 'X', 'N', 'L', 'S']

# Feedback analysis
# Guess 1: 95WC
# - One number correct but wrong position: 9 or 5
# - One letter correct and correct position: W
# - One letter incorrect and too early: C

# Guess 2: 12XN
# - One number correct and correct position: 2
# - Both letters incorrect: X, N

# Guess 3: 92LS
# - One number correct but wrong position: 9 or 2
# - One letter correct but wrong position: L

# Deduce the correct combination
correct_numbers = ['2', '9']
correct_letters = ['W', 'L']

# Verify the deduction
print(correct_numbers + correct_letters)