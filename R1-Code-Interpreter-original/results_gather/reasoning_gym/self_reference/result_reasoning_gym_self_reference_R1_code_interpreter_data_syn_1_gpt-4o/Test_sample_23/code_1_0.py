def check_statements():
    possible_solutions = 0
    for T in range(8):
        F = 7 - T
        # Check each statement
        s1 = (T >= 3)
        s2 = (F <= 1)
        s3 = (T == 6)
        s4 = (F == 2)
        s5 = (s3 != s4)
        s6 = (T in [2, 3, 5, 7])
        s7 = (F in [4, 6, 8, 9, 10, 12, 14])  # Composite numbers

        # Check if all statements are consistent
        if s1 and s2 and s3 and s4 and s5 and s6 and s7:
            possible_solutions += 1

    print(possible_solutions)

check_statements()