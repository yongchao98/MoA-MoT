def check_guess(numbers, letters):
    # Check first guess (97EJ)
    if numbers[0] >= 9 or numbers[1] >= 7:
        return False
    
    # Check second guess (43OM)
    if numbers[0] >= 4 or numbers[1] >= 3:
        return False
    
    # Check third guess (20BQ)
    num_correct = 0
    if numbers[0] == 2 or numbers[0] == 0:
        num_correct += 1
    if numbers[1] == 2 or numbers[1] == 0:
        num_correct += 1
    if num_correct != 1:
        return False
    if max(numbers) >= 2 and max(numbers) != 2 and max(numbers) != 0:
        return False
    
    # Letter checks
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # From first guess (97EJ)
    ej_correct = 0
    if 'E' in letters or 'J' in letters:
        ej_correct += 1
    if ej_correct != 1:
        return False
    
    # From second guess (43OM)
    if any(l >= 'O' for l in letters):
        return False
    
    # From third guess (20BQ)
    bq_correct = 0
    if 'B' in letters or 'Q' in letters:
        bq_correct += 1
    if bq_correct != 1:
        return False
    if 'Q' in letters:  # Q is too late in alphabet
        return False
    
    return True

# Generate all possible combinations
valid_solutions = []
numbers = list(range(10))
letters = [chr(i) for i in range(65, 91)]  # A-Z

for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        if check_guess([n1, n2], [l1, l2]):
                            valid_solutions.append([str(n1), str(n2), l1, l2])

print(valid_solutions)