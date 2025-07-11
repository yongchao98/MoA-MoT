def solve_statements():
    """
    Solves the graph theory statements by applying a derived condition.

    The analysis shows that a statement of the form A = O(B/d^x + nd^y) is true
    if and only if the exponents x and y satisfy the condition y >= x + 1.
    """

    statements = [
        # Statement: A = O(B/d^x + nd^y)
        # Parameters: (statement_number, x, y)
        (1, 0.5, 2.0),
        (2, 0.99, 2.0),
        (3, 1.0, 2.0),
        (4, 1.5, 2.0),
        (5, 1.3, 2.6),
        (6, 1.0, 1.5),
    ]

    final_answer = ""
    print("Checking each statement based on the condition y >= x + 1:")
    
    for num, x, y in statements:
        is_true = (y >= x + 1)
        result_char = 'Y' if is_true else 'N'
        final_answer += result_char
        
        # Output the numbers and logic for each statement
        print(f"Statement {num}: A = O(B/d^{x} + nd^{y})")
        print(f"  x={x}, y={y}")
        print(f"  Checking if {y} >= {x} + 1: {y} >= {x+1}")
        print(f"  Result: {'True' if is_true else 'False'}")
        print("-" * 20)
        
    print(f"Final Answer String: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_statements()
