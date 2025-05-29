from itertools import permutations

def check_guess(guess, actual):
    # Check numbers
    num_correct_pos = sum(1 for i in range(2) if guess[i] == actual[i])
    num_correct_wrong_pos = sum(1 for i in range(2) for j in range(2) if i != j and guess[i] == actual[j])
    
    # Check letters
    let_correct_pos = sum(1 for i in range(2,4) if guess[i] == actual[i])
    let_too_late = sum(1 for i in range(2,4) if guess[i] > actual[i])
    
    return num_correct_pos, num_correct_wrong_pos, let_correct_pos, let_too_late

def matches_condition(candidate):
    # Convert candidate to string format for comparison
    cand_str = ''.join(map(str, candidate))
    
    # Check all conditions
    conditions = [
        ('19HE', (0, 0, 1, 1)),
        ('75IS', (0, 0, 0, 2)),
        ('58JW', (0, 1, 0, 2)),
        ('64XZ', (0, 0, 0, 2)),
        ('07HM', (0, 1, 0, 2)),
        ('35WS', (0, 0, 0, 2)),
        ('47HL', (0, 0, 0, 2)),
        ('23QE', (0, 0, 1, 1)),
        ('48AJ', (0, 1, 1, 1))
    ]
    
    for guess, expected in conditions:
        result = check_guess(guess, cand_str)
        if result != expected:
            return False
    return True

# Generate all possible combinations
digits = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

all_solutions = []
for d1 in digits:
    for d2 in digits:
        if d1 != d2:  # no repeating numbers
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # no repeating letters
                        candidate = [d1, d2, l1, l2]
                        if matches_condition(candidate):
                            all_solutions.append(candidate)

print(all_solutions)