def is_too_large(guess, actual):
    return int(guess) > int(actual)

def is_too_early(guess, actual):
    return ord(guess) < ord(actual)

def check_condition1(solution):
    # 68BC: both numbers too large, both letters too early
    return (is_too_large('6', solution[0]) and 
            is_too_large('8', solution[1]) and 
            is_too_early('B', solution[2]) and 
            is_too_early('C', solution[3]))

def check_condition2(solution):
    # 41WA: one number correct wrong pos, one too large, both letters incorrect
    nums = [solution[0], solution[1]]
    return (('1' in nums and '1' != solution[0]) or ('4' in nums and '4' != solution[1])) and \
           (is_too_large('4', solution[0]) or is_too_large('1', solution[1])) and \
           solution[2] != 'W' and solution[3] != 'A'

def check_condition3(solution):
    # 42MQ: both numbers too large, one letter wrong pos, one incorrect
    letters = [solution[2], solution[3]]
    return (is_too_large('4', solution[0]) and 
            is_too_large('2', solution[1]) and 
            ('M' in letters or 'Q' in letters) and 
            'M' != solution[2] and 'Q' != solution[3])

def check_condition4(solution):
    # 59CR: both numbers too large, both letters incorrect
    return (is_too_large('5', solution[0]) and 
            is_too_large('9', solution[1]) and 
            solution[2] != 'C' and solution[3] != 'R')

def check_condition5(solution):
    # 12AE: one number correct position, one too large, one letter correct position, one too early
    return ((solution[0] == '1' or solution[1] == '2') and 
            (is_too_large('2', solution[0]) or is_too_large('1', solution[1])) and
            (solution[2] == 'A' or solution[3] == 'E') and 
            (is_too_early('E', solution[2]) or is_too_early('A', solution[3])))

# Generate and test all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

found = False
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        solution = [n1, n2, l1, l2]
                        if (check_condition1(solution) and 
                            check_condition2(solution) and 
                            check_condition3(solution) and 
                            check_condition4(solution) and 
                            check_condition5(solution)):
                            print(solution)
                            found = True
                            break
                if found:
                    break
            if found:
                break
        if found:
            break