# Define the statements
def statement_1(true_count):
    return true_count >= 7

def statement_2(false_count):
    return false_count <= 3

def statement_3(true_count):
    return true_count == 4

def statement_4(false_count):
    return false_count == 3

def statement_5(true_count):
    return (statement_3(true_count) and not statement_4(7 - true_count)) or (not statement_3(true_count) and statement_4(7 - true_count))

def statement_6(true_count):
    return true_count in [2, 3, 5, 7]

def statement_7(false_count):
    return false_count in [4, 6]

# Check all possible true statement counts
possible_solutions = 0
for true_count in range(8):
    false_count = 7 - true_count
    if (not statement_1(true_count) and
        statement_2(false_count) and
        statement_5(true_count) and
        statement_6(true_count) and
        statement_7(false_count)):
        possible_solutions += 1

print(possible_solutions)