def check_guess(numbers, letters):
    # Check first guess (97EJ)
    if numbers[0] >= 9 or numbers[1] >= 7:
        return False
    
    # Check second guess (43OM)
    if numbers[0] >= 4 or numbers[1] >= 3:
        return False
    
    # Check third guess (20BQ)
    # 0 must be in the solution, 2 is too large
    if 2 in numbers or 0 not in numbers:
        return False
    
    # Letter checks
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # From first guess (97EJ)
    # Since B must be in solution (from third guess analysis),
    # E or J must be in the correct position
    if not (('E' in letters and letters.index('E') == letters.index('E')) or 
            ('J' in letters and letters.index('J') == letters.index('J'))):
        return False
    
    # From second guess (43OM)
    if any(l >= 'O' for l in letters):
        return False
    
    # From third guess (20BQ)
    if 'B' not in letters:
        return False
    if 'Q' in letters:  # Q is too late in alphabet
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
numbers = list(range(10))
letters = [chr(i) for i in range(65, 91)]  # A-Z

for n1 in range(3):  # Numbers must be less than 3
    for n2 in range(3):
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        if check_guess([n1, n2], [l1, l2]):
                            valid_solutions.append([str(n1), str(n2), l1, l2])

print(valid_solutions)