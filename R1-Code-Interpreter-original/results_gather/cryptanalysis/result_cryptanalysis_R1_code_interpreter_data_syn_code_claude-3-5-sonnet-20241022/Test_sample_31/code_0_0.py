from itertools import permutations

def check_guess(guess):
    # All conditions from the hints
    
    # Condition 1: 92BJ - both numbers wrong, both letters wrong
    if guess[0] in ['9', '2'] or guess[1] in ['9', '2'] or guess[2] in ['B', 'J'] or guess[3] in ['B', 'J']:
        return False
    
    # Condition 2: 06RU - one number correct in position, one too small, both letters too late
    if not ((guess[0] == '0' and int(guess[1]) > 6) or (guess[0] == '6' and int(guess[1]) > 0)):
        return False
    if guess[2] >= 'R' or guess[3] >= 'R':
        return False
    
    # Condition 3: 84CV - both numbers wrong, one letter correct in position, one too late
    if guess[0] in ['8', '4'] or guess[1] in ['8', '4']:
        return False
    if not ((guess[2] == 'C' and guess[3] < 'V') or (guess[2] < 'V' and guess[3] == 'C')):
        return False
    
    # Condition 4: 75VY - one number correct in position, one too small
    if not ((guess[0] == '7' and int(guess[1]) > 5) or (guess[0] == '5' and int(guess[1]) > 7)):
        return False
    
    # Condition 5: 16UD - one number correct in position, one too small
    if not ((guess[0] == '1' and int(guess[1]) > 6) or (guess[0] == '6' and int(guess[1]) > 1)):
        return False
    if guess[2] in ['U', 'D'] or guess[3] in ['U', 'D']:
        return False
    
    # Condition 6: 98NG - both numbers too large
    if int(guess[0]) >= 9 or int(guess[1]) >= 8:
        return False
    if guess[2] in ['N', 'G'] or guess[3] in ['N', 'G']:
        return False
    
    # Condition 7: 03QD - both numbers too small, one letter correct wrong position
    if int(guess[0]) <= 0 or int(guess[1]) <= 3:
        return False
    if not ((guess[2] == 'Q' and guess[3] != 'D') or (guess[3] == 'Q' and guess[2] != 'D')):
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

valid_solutions = []
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        guess = [n1, n2, l1, l2]
                        if check_guess(guess):
                            valid_solutions.append(guess)

print(valid_solutions)