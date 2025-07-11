def solve_profit_expression():
    """
    This function determines the algebraic expression for a firm's profit
    based on the areas S, T, and H as defined in the problem.
    """

    # 1. The fundamental formula for profit is Total Revenue (TR) minus Total Cost (TC).
    # Profit = TR - TC

    # 2. Identify Total Revenue (TR).
    # TR is Price (P) times quantity (q1).
    # The area S is defined as a rectangle with area P * q1.
    # Therefore, TR is represented by the area S.
    tr_expression = "S"

    # 3. Identify Total Cost (TC).
    # TC is Average Total Cost (ATC) times quantity (q1).
    # The area H is defined as a rectangle with area ATC * q1.
    # Therefore, TC is represented by the area H.
    tc_expression = "H"

    # 4. The final expression for profit is TR - TC.
    profit_expression_parts = [tr_expression, "-", tc_expression]

    # Print the step-by-step derivation
    print("Step 1: The firm's profit is defined as Total Revenue (TR) - Total Cost (TC).")
    print(f"Step 2: Total Revenue (TR) is represented by the area {tr_expression}.")
    print(f"Step 3: Total Cost (TC) is represented by the area {tc_expression}.")
    print("\nStep 4: Substituting these into the profit formula gives the final expression.")
    print("The area representing the firm's profit is:")
    
    # Print the final equation with each component separated as requested.
    print(profit_expression_parts[0], profit_expression_parts[1], profit_expression_parts[2])

solve_profit_expression()