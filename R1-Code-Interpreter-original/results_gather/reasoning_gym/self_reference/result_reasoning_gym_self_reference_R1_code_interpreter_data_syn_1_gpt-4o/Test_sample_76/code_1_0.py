from sympy import isprime

def check_statements():
    possible_solutions = 0
    for true_count in range(2, 8):  # At least 2 true statements
        false_count = 7 - true_count
        if true_count == 5 and isprime(true_count) and false_count in [4, 6, 8]:
            possible_solutions += 1
    print(possible_solutions)

check_statements()