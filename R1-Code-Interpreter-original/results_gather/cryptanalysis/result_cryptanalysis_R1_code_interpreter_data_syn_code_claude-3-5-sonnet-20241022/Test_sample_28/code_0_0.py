def check_conditions(guess):
    # All numbers that can't be in the solution
    forbidden_nums = {1,5,8,0,3,4}
    # All letters that can't be in the solution
    forbidden_letters = {'Z','Q','I','K','U','M','J','Y','D'}
    # Letters that are too early in alphabet
    early_letters = {'I','K'}
    
    n1, n2 = int(guess[0]), int(guess[1])
    l1, l2 = guess[2], guess[3]
    
    # Check forbidden numbers
    if n1 in forbidden_nums or n2 in forbidden_nums:
        return False
    
    # Check forbidden letters
    if l1 in forbidden_letters or l2 in forbidden_letters:
        return False
    
    # Check condition from guess 3
    if n1 != 2 and n2 != 2:
        return False
        
    # Check condition from guess 6
    if n1 != 6:  # 6 must be in first position
        return False
    if 'R' not in [l1, l2]:  # R must be present
        return False
        
    # Generate all possible solutions
    return True

# Generate all possible solutions
solutions = []
numbers = range(10)
letters = [chr(i) for i in range(65, 91)]  # A-Z

for n1 in numbers:
    for n2 in numbers:
        if n1 == n2:  # numbers can't repeat
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:  # letters can't repeat
                    continue
                guess = f"{n1}{n2}{l1}{l2}"
                if check_conditions(guess):
                    solutions.append([str(n1), str(n2), l1, l2])

print(solutions)