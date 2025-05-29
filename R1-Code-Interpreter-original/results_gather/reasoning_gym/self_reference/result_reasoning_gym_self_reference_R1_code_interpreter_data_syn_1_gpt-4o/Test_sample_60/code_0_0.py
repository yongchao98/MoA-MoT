def check_statements():
    possible_solutions = 0
    for T in range(8):  # T can be from 0 to 7
        F = 7 - T
        s1 = (T >= 4)
        s2 = (F <= 0)
        s3 = (T == 4)
        s4 = (F == 3)
        s5 = (s3 != s4)
        s6 = (T in [2, 3, 5, 7])
        s7 = (F in [4, 6])
        
        if all([s1, s2, s5, s6]) and not s7:
            possible_solutions += 1
    
    print(possible_solutions)

check_statements()