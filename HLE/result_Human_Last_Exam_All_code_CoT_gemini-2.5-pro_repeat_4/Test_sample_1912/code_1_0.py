def solve_sequence():
    """
    Solves the sequence puzzle by identifying and applying a two-part pattern.
    The pattern for 'a + b' depends on the ratio of b to a.
    """
    equations = [(1, 2), (2, 4), (3, 9), (5, 10), (6, 12)]
    final_answer = 0

    print("Analyzing the sequence with the discovered pattern:\n")

    for a, b in equations:
        # Check if b is 2 times a
        if b == 2 * a:
            # The pattern is: a * b + (a - 2)
            multiplication_part = a * b
            adjustment_part = a - 2
            result = multiplication_part + adjustment_part
            
            # Print the equation with its components
            if adjustment_part < 0:
                print(f"{a} + {b} = {multiplication_part} - {abs(adjustment_part)} = {result}")
            else:
                print(f"{a} + {b} = {multiplication_part} + {adjustment_part} = {result}")
            
            final_answer = result
        
        # Check if b is 3 times a
        elif b == 3 * a:
            # The pattern is simply: a * b
            result = a * b
            print(f"{a} + {b} = {a} * {b} = {result}")
            final_answer = result

# Execute the function to find and print the solutions
solve_sequence()