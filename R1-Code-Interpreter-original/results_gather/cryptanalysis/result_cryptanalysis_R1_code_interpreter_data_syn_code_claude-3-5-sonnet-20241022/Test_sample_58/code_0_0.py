from itertools import permutations

def check_guess(guess):
    # All known conditions
    conditions = [
        # 10PR - both numbers wrong and too small
        lambda g: not (g[0] == '1' or g[1] == '0') and int(g[0]) > 1 and int(g[1]) > 0,
        # 79EN - one number correct but wrong position, one too large
        lambda g: (g[0] == '9' or g[1] == '9' or g[0] == '7' or g[1] == '7') and not (g[0] == '7' and g[1] == '9') and not (g[0] == '9' and g[1] == '7'),
        # 70GE - one number correct but wrong position, one too small
        lambda g: (g[0] == '7' or g[1] == '7' or g[0] == '0' or g[1] == '0') and not (g[0] == '7' and g[1] == '0') and not (g[0] == '0' and g[1] == '7'),
        # 94QE - both numbers wrong
        lambda g: g[0] != '9' and g[1] != '4',
        # 05LE - one number correct but wrong position, one too small
        lambda g: (g[0] == '0' or g[1] == '0' or g[0] == '5' or g[1] == '5') and not (g[0] == '0' and g[1] == '5') and not (g[0] == '5' and g[1] == '0'),
        # 83CS - both numbers wrong
        lambda g: g[0] != '8' and g[1] != '3',
        # 03ZF - both numbers too small
        lambda g: int(g[0]) > 0 and int(g[1]) > 3
    ]
    
    # Check all number conditions
    return all(cond(guess) for cond in conditions)

def check_letters(letters):
    # Letter conditions
    if 'Z' not in letters[1]:  # Z must be in second position
        return False
    # Letters that can't be used
    invalid_letters = {'P', 'R', 'E', 'N', 'G', 'Q', 'L', 'C', 'S', 'F'}
    if any(l in invalid_letters for l in letters):
        return False
    return True

# Generate all possible number combinations
numbers = [str(i) for i in range(10)]
valid_numbers = []
for num1, num2 in permutations(numbers, 2):
    if check_guess([num1, num2]):
        valid_numbers.append([num1, num2])

# Generate all possible letter combinations
letters = [chr(i) for i in range(65, 91)]  # A-Z
valid_combinations = []
for num in valid_numbers:
    for l1 in letters:
        for l2 in letters:
            if l1 != l2 and check_letters([l1, l2]):
                valid_combinations.append(num + [l1, 'Z'])

if valid_combinations:
    print(valid_combinations)