def check_conditions(guess):
    numbers = guess[:2]
    letters = guess[2:]
    
    # Check number conditions
    if '2' in numbers or '0' in numbers: return False  # from 20IX
    if '7' in numbers or '8' in numbers: return False  # from 78LK
    if '4' in numbers or '8' in numbers: return False  # from 48TP
    if '2' in numbers or '5' in numbers: return False  # from 25LH
    if '7' in numbers or '2' in numbers: return False  # from 72NO
    if '6' in numbers or '5' in numbers: return False  # from 65RX
    
    # First digit must be 1 (from 21XS)
    if numbers[0] != '1': return False
    
    # Check letter conditions
    if 'O' not in letters[1]: return False  # O must be second letter
    
    # Letters must be between L-H and T-X in alphabet
    for letter in letters:
        if letter < 'H' or letter > 'T': return False
        if letter in 'IXYTLKLDCOXSTPLHNORX': return False
    
    return True

# Generate all possible combinations
valid_combinations = []
for n1 in '1':  # first number must be 1
    for n2 in '0123456789':
        if n1 != n2:  # numbers can't repeat
            for l1 in 'IJKLMNOPQRST':
                for l2 in 'O':  # second letter must be O
                    guess = n1 + n2 + l1 + l2
                    if check_conditions(guess):
                        valid_combinations.append([n1, n2, l1, l2])

print(valid_combinations)