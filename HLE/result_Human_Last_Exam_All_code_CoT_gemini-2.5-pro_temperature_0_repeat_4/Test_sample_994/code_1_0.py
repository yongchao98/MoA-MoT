import itertools
import operator

def solve_24_game():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can form the number 24.
    """

    def can_make_24(nums):
        """
        Checks if a given list of 4 numbers can result in 24 using +, -, *, /.
        """
        # Use a small tolerance for floating point comparisons
        epsilon = 1e-6
        
        # Try all unique permutations of the numbers
        for p_nums in set(itertools.permutations(nums)):
            a, b, c, d = p_nums
            
            # Try all 4^3 = 64 combinations of operators
            for op1, op2, op3 in itertools.product([operator.add, operator.sub, operator.mul, operator.truediv], repeat=3):
                
                # Parenthesization 1: ((a op1 b) op2 c) op3 d
                try:
                    res1 = op3(op2(op1(a, b), c), d)
                    if abs(res1 - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass

                # Parenthesization 2: (a op1 (b op2 c)) op3 d
                try:
                    res2 = op3(op1(a, op2(b, c)), d)
                    if abs(res2 - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass

                # Parenthesization 3: a op1 ((b op2 c) op3 d)
                try:
                    res3 = op1(a, op3(op2(b, c), d))
                    if abs(res3 - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass

                # Parenthesization 4: a op1 (b op2 (c op3 d))
                try:
                    res4 = op1(a, op2(b, op3(c, d)))
                    if abs(res4 - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass

                # Parenthesization 5: (a op1 b) op2 (c op3 d)
                try:
                    # Check for division by zero in the second part
                    val2 = op3(c, d)
                    res5 = op2(op1(a, b), val2)
                    if abs(res5 - 24) < epsilon:
                        return True
                except ZeroDivisionError:
                    pass
        
        return False

    # 1. Generate all possible combinations of four cards (values 1-10)
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    # 2. Count how many combinations can make 24
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # 3. Calculate and print the final percentage
    percentage = solvable_count / total_combinations
    
    # Print the final equation as requested
    print(f"{solvable_count} / {total_combinations} = {percentage:.4f}")

solve_24_game()
<<<0.7720>>>