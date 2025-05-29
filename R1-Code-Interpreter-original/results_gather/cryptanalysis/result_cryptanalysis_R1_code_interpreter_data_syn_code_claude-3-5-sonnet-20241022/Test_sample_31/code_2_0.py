def verify_solution(guess):
    n1, n2, l1, l2 = guess
    
    # Additional verification for Condition 2: 06RU
    # We know one number is correct and in correct position
    # Since both solutions have 7,6, let's verify which position matches with 06
    if n1 != '0' and n2 != '6':  # neither matches position with 06
        return False
    
    # Additional verification for Condition 4: 75VY
    # We know one number is correct and in correct position
    # This means 7 must be in position 1 (because if 5 was correct, it would match with solution)
    if n1 != '7':
        return False
    
    # Additional verification for Condition 3: 84CV
    # We know one letter is correct and in correct position
    # Since C and Q are our letters, C must be the one in correct position
    # V is too late in alphabet
    if l1 != 'C':  # C must be in first letter position
        return False
    
    return True

solutions = [['7', '6', 'C', 'Q'], ['7', '6', 'Q', 'C']]
final_solutions = [sol for sol in solutions if verify_solution(sol)]
print(final_solutions)