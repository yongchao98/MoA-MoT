def check_guess(guess, target):
    # Returns feedback similar to the given conditions
    num_pos = 0  # numbers correct and in position
    num_wrong = 0  # numbers incorrect
    num_wrong_larger = 0  # numbers too large
    let_pos = 0  # letters correct and in position
    let_wrong = 0  # letters correct but wrong position
    
    # Check numbers
    for i in range(2):
        if guess[i] == target[i]:
            num_pos += 1
        elif guess[i] > target[i]:
            num_wrong_larger += 1
        else:
            num_wrong += 1
            
    # Check letters
    for i in range(2, 4):
        if guess[i] == target[i]:
            let_pos += 1
        elif guess[i] in target[2:]:
            let_wrong += 1
    
    return (num_pos, num_wrong_larger, num_wrong, let_pos, let_wrong)

def matches_condition(guess, target, condition):
    feedback = check_guess(guess, target)
    return feedback == condition

def find_password():
    numbers = range(10)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # Store all conditions
    conditions = [
        ([8,0,'Z','V'], (0,0,2,0,0)),  # both numbers wrong, both letters wrong
        ([9,1,'P','E'], (1,1,0,0,0)),  # one number correct position, one too large
        ([5,1,'D','H'], (1,0,1,0,0)),  # one number correct position
        ([1,9,'H','C'], (0,1,1,1,0)),  # one number wrong pos, one too large, one letter correct pos
        ([3,8,'C','A'], (0,0,2,0,1)),  # both numbers wrong, one letter wrong pos
        ([2,5,'F','G'], (0,0,2,0,0)),  # both wrong
        ([2,1,'G','K'], (1,0,1,0,0)),  # one number correct pos
        ([4,8,'F','L'], (0,0,2,0,0)),  # both wrong
        ([8,0,'D','L'], (0,0,2,0,0)),  # both wrong
        ([4,3,'Z','O'], (0,0,2,0,0)),  # both wrong
        ([9,2,'S','C'], (0,0,2,1,0)),  # both numbers wrong, one letter correct pos
        ([0,3,'X','O'], (0,0,2,1,0)),  # both numbers wrong, one letter correct pos
        ([2,3,'J','H'], (0,0,2,0,0)),  # both wrong
        ([0,5,'T','Q'], (0,0,2,0,0)),  # both wrong
        ([6,0,'Q','F'], (0,0,2,0,0))   # both wrong
    ]
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    
                    target = [n1, n2, l1, l2]
                    valid = True
                    
                    # Check if this combination satisfies all conditions
                    for guess, condition in conditions:
                        if not matches_condition(guess, target, condition):
                            valid = False
                            break
                    
                    if valid:
                        print(f"Found possible solution: {target}")

find_password()