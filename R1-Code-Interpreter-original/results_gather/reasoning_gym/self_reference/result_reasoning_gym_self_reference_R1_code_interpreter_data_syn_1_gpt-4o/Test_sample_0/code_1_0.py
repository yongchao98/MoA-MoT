def check_statements():
    possible_solutions = 0
    for T in range(8):  # T can be from 0 to 7
        F = 7 - T
        # Check each statement
        s1 = (T >= 3)
        s2 = (F <= 1)
        s3 = (T == 3)
        s4 = (F == 5)
        s5 = (s3 != s4)
        s6 = (T in [2, 3, 5, 7])  # Prime numbers
        s7 = (F in [4, 6, 8])  # Composite numbers

        # Count the number of true statements
        true_statements = sum([s1, s2, s3, s4, s5, s6, s7])

        # Check if the number of true statements matches T
        if true_statements == T:
            possible_solutions += 1

    print(possible_solutions)

check_statements()