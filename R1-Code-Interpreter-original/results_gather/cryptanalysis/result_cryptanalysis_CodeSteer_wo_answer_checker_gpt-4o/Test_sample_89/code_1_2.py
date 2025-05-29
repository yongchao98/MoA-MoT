from itertools import permutations

def is_valid_combination(combination):
    num1, num2, let1, let2 = combination
    
    # Feedback conditions
    conditions = [
        # Guess: 16TE
        (num1 not in [1, 6] and num2 not in [1, 6] and let1 not in ['T', 'E'] and let2 not in ['T', 'E']),
        # Guess: 72QY
        (num1 == 7 and num2 > 2 and ((let1 == 'Q' and let2 != 'Y') or (let1 != 'Q' and let2 == 'Q'))),
        # Guess: 23AX
        (num1 == 3 and num2 > 2 and let1 not in ['A', 'X'] and let2 not in ['A', 'X']),
        # Guess: 71KQ
        (num1 == 7 and num2 > 1 and let1 not in ['K', 'Q'] and let2 not in ['K', 'Q']),
        # Guess: 91TW
        (num1 not in [9, 1] and num2 not in [9, 1] and let1 not in ['T', 'W'] and let2 not in ['T', 'W']),
        # Guess: 51PL
        (num1 not in [5, 1] and num2 not in [5, 1] and ((let1 == 'L' and let2 != 'P') or (let1 != 'P' and let2 == 'L')))
    ]
    
    return all(conditions)

def find_password():
    numbers = range(10)
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for num_comb in permutations(numbers, 2):
        for let_comb in permutations(letters, 2):
            combination = (*num_comb, *let_comb)
            if is_valid_combination(combination):
                print(f"Valid combination found: {combination}")
                return combination

    print("No valid combination found.")
    return None

password = find_password()
if password:
    print(f"<<< {list(map(str, password))} >>>")
else:
    print("Failed to find a valid password.")