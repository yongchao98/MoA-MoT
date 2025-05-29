# Define the statements as functions
def statement_1(truths):
    return truths >= 7

def statement_2(truths):
    return truths >= 5

def statement_3(truths):
    return truths == 1

def statement_4(truths):
    return truths == 2

def statement_5(truths):
    return (statement_3(truths) or statement_4(truths)) and not (statement_3(truths) and statement_4(truths))

def statement_6(truths):
    return truths in {2, 3, 5, 7}  # Prime numbers

def statement_7(falses):
    return falses in {4, 6}  # Composite numbers

# Check all combinations of true statements
possible_solutions = 0
for truths in range(8):
    falses = 7 - truths
    if (not statement_1(truths) and
        statement_2(truths) and
        not statement_3(truths) and
        statement_4(truths) and
        statement_5(truths) and
        statement_6(truths) and
        statement_7(falses)):
        possible_solutions += 1

print(possible_solutions)