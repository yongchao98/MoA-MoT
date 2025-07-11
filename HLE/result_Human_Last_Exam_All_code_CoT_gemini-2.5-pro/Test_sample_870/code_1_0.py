import itertools

def solve_24_game(numbers):
    """
    Finds and prints a solution to the 24-point game for a given list of four numbers.
    """
    # Using a set gets unique permutations, which is efficient for repeated numbers.
    num_permutations = set(itertools.permutations(numbers))
    operators = ['+', '-', '*', '/']
    op_combinations = list(itertools.product(operators, repeat=3))

    # Pattern 1: ((n1 op1 n2) op2 n3) op3 n4
    for nums in num_permutations:
        n1, n2, n3, n4 = nums
        for op1, op2, op3 in op_combinations:
            try:
                # Use floating point numbers for calculation to handle division properly
                val1 = eval(f"{n1} {op1} {n2}")
                val2 = eval(f"{val1} {op2} {n3}")
                final_val = eval(f"{val2} {op3} {n4}")

                # Check if the result is 24 (with a small tolerance for float inaccuracies)
                if abs(final_val - 24) < 0.0001:
                    print("A solution has been found! Here is the calculation:")
                    print(f"Step 1: {n1} {op1} {n2} = {val1}")
                    print(f"Step 2: {val1} {op2} {n3} = {val2}")
                    print(f"Step 3: {val2} {op3} {n4} = {int(final_val)}")
                    print("\nThe final equation is:")
                    print(f"(({n1} {op1} {n2}) {op2} {n3}) {op3} {n4} = 24")
                    return
            except ZeroDivisionError:
                continue

    # Pattern 2: (n1 op1 n2) op2 (n3 op3 n4)
    for nums in num_permutations:
        n1, n2, n3, n4 = nums
        for op1, op2, op3 in op_combinations:
            try:
                # Avoid division by zero in the second term
                if op3 == '/' and n4 == 0:
                    continue
                
                val1 = eval(f"{n1} {op1} {n2}")
                val2 = eval(f"{n3} {op3} {n4}")
                final_val = eval(f"{val1} {op2} {val2}")

                if abs(final_val - 24) < 0.0001:
                    print("A solution has been found! Here is the calculation:")
                    print(f"Step 1: {n1} {op1} {n2} = {val1}")
                    print(f"Step 2: {n3} {op3} {n4} = {val2}")
                    print(f"Step 3: {val1} {op2} {val2} = {int(final_val)}")
                    print("\nThe final equation is:")
                    print(f"({n1} {op1} {n2}) {op2} ({n3} {op3} {n4}) = 24")
                    return
            except ZeroDivisionError:
                continue
    
    print("No solution was found for the given numbers.")

# The numbers for this puzzle
puzzle_numbers = [4, 4, 10, 10]
solve_24_game(puzzle_numbers)