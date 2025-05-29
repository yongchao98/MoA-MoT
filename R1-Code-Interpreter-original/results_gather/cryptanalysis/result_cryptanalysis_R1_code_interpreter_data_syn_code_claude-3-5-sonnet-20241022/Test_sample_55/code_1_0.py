def check_additional_constraints(solutions):
    valid = []
    for sol in solutions:
        first_num = int(sol[0])
        valid_solution = True
        
        # From guess 1 (43AJ): 3 must be too small compared to first number
        if first_num <= 3:
            valid_solution = False
            
        # From guess 2 (40HZ): 0 must be too small compared to first number
        if first_num <= 0:
            valid_solution = False
            
        # From guess 7 (29FP): 2 and 9 must be wrong
        if first_num in [2, 9]:
            valid_solution = False
            
        if valid_solution:
            valid.append(sol)
    
    return valid

initial_solutions = [['5', '4', 'Y', 'Q'], ['6', '4', 'Y', 'Q'], ['7', '4', 'Y', 'Q'], ['9', '4', 'Y', 'Q']]
final_solutions = check_additional_constraints(initial_solutions)
print(final_solutions)