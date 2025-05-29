from itertools import product
import string

def check_guess(guess, actual):
    # Returns feedback similar to the given conditions
    numbers_correct = [0, 0]  # [correct position, wrong position]
    letters_correct = [0, 0]  # [correct position, wrong position]
    
    for i in range(2):
        if guess[i] == actual[i]:
            numbers_correct[0] += 1
        elif guess[i] in actual[:2]:
            numbers_correct[1] += 1
            
    for i in range(2,4):
        if guess[i] == actual[i]:
            letters_correct[0] += 1
        elif guess[i] in actual[2:]:
            letters_correct[1] += 1
            
    return numbers_correct, letters_correct

def is_consistent(candidate):
    # Check all conditions
    conditions = [
        ('79AG', ([0,1], [0,0])),
        ('15BN', ([1,0], [0,0])),
        ('46QN', ([0,0], [0,0])),
        ('16BC', ([0,0], [0,0])),
        ('81BA', ([0,0], [0,0])),
        ('69IH', ([0,1], [0,0])),
        ('05WO', ([1,0], [0,0])),
        ('74PV', ([0,0], [1,0])),
        ('26BF', ([0,0], [0,0])),
        ('63TH', ([0,0], [0,0])),
        ('74HT', ([0,0], [0,0])),
        ('63OY', ([0,0], [0,0])),
        ('06XR', ([0,0], [0,0])),
        ('63HS', ([0,0], [0,0])),
        ('13EX', ([0,0], [1,0]))
    ]
    
    for guess, expected in conditions:
        result = check_guess(guess, candidate)
        if result != expected:
            return False
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = list(string.ascii_uppercase)

valid_candidates = []
for n1, n2 in product(numbers, numbers):
    if n1 != n2:  # no repeating numbers
        for l1, l2 in product(letters, letters):
            if l1 != l2:  # no repeating letters
                candidate = n1 + n2 + l1 + l2
                if is_consistent(candidate):
                    valid_candidates.append(candidate)

print(valid_candidates)