from itertools import combinations

# Define the statements as functions
def statement_1(true_count):
    return true_count >= 4

def statement_2(false_count):
    return false_count <= 1

def statement_3(true_count):
    return true_count == 4

def statement_4(false_count):
    return false_count == 3

def statement_5(s3, s4):
    return (s3 or s4) and not (s3 and s4)

def statement_6(true_count):
    return true_count in {2, 3, 5, 7}

def statement_7(false_count):
    return false_count in {4, 6, 8, 9, 10}

# Check all possible combinations of true/false for 7 statements
possible_solutions = 0
for true_statements in combinations(range(7), 3):
    true_count = len(true_statements)
    false_count = 7 - true_count
    
    s1 = statement_1(true_count)
    s2 = statement_2(false_count)
    s3 = statement_3(true_count)
    s4 = statement_4(false_count)
    s5 = statement_5(s3, s4)
    s6 = statement_6(true_count)
    s7 = statement_7(false_count)
    
    # Check if the configuration is consistent
    if sum([s1, s2, s3, s4, s5, s6, s7]) == 3:
        possible_solutions += 1

print(possible_solutions)