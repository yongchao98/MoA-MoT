def check_guess(guess, actual):
    # Returns feedback based on the rules
    correct_nums = 0
    wrong_pos_nums = 0
    correct_letters = 0
    wrong_pos_letters = 0
    
    # Check numbers (first two positions)
    for i in range(2):
        if guess[i] == actual[i]:
            correct_nums += 1
        elif guess[i] in actual[:2]:
            wrong_pos_nums += 1
            
    # Check letters (last two positions)
    for i in range(2,4):
        if guess[i] == actual[i]:
            correct_letters += 1
        elif guess[i] in actual[2:]:
            wrong_pos_letters += 1
            
    return (correct_nums, wrong_pos_nums, correct_letters, wrong_pos_letters)

def is_valid_solution(solution):
    # Check all given conditions
    if check_guess(['6','8','B','C'], solution) != (0,0,0,0):  # both numbers too large, letters too early
        return False
    if check_guess(['4','1','W','A'], solution) != (0,1,0,0):  # one number wrong pos, one too large
        return False
    if check_guess(['4','2','M','Q'], solution) != (0,0,0,1):  # numbers too large, one letter wrong pos
        return False
    if check_guess(['5','9','C','R'], solution) != (0,0,0,0):  # all incorrect
        return False
    if check_guess(['1','2','A','E'], solution) != (1,0,1,0):  # one number correct pos, one letter correct pos
        return False
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

best_solution = None

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        solution = [n1, n2, l1, l2]
                        if is_valid_solution(solution):
                            best_solution = solution
                            break
                if best_solution:
                    break
            if best_solution:
                break
        if best_solution:
            break

print(best_solution)