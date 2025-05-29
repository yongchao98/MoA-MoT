from sympy import isprime

def check_statements():
    possible_solutions = 0
    for true_count in range(1, 8):
        false_count = 7 - true_count
        # Check conditions
        if (true_count >= 4 and
            (true_count == 4 or true_count == 3) and
            (true_count == 4) != (false_count == 3) and
            isprime(true_count) and
            false_count in [4, 6]):
            possible_solutions += 1
    print(possible_solutions)

check_statements()