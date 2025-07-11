def solve_zoo_puzzle():
    """
    This function solves the "Fun Facts From The Zoo" puzzle by printing
    the hidden equation in each line.
    """
    equations = [
        {"num1": 10, "op": "-", "num2": 5, "result": 5},
        {"num1": 2, "op": "+", "num2": 3, "result": 5},
        {"num1": 2, "op": "**", "num2": 3, "result": 8},
        {"num1": 20, "op": "/", "num2": 4, "result": 5},
        {"num1": 1, "op": "*", "num2": 7, "result": 7},
        {"num1": 1, "op": "+", "num2": 4, "result": 5},
        {"num1": 6, "op": "/", "num2": 2, "result": 3},
        {"num1": 9, "op": "-", "num2": 6, "result": 3},
        {"num1": 12, "op": "/", "num2": 3, "result": 4},
        {"num1": 8, "op": "-", "num2": 2, "result": 6},
        {"num1": 9, "op": "-", "num2": 2, "result": 7},
        {"num1": 4, "op": "*", "num2": 2, "result": 8},
        {"num1": 8, "op": "-", "num2": 3, "result": 5},
        {"num1": 0, "op": "+", "num2": 4, "result": 4},
        {"num1": 8, "op": "-", "num2": 2, "result": 6},
        {"num1": 0, "op": "+", "num2": 4, "result": 4},
        {"num1": 9, "op": "-", "num2": 4, "result": 5},
        {"num1": 16, "op": "/", "num2": 2, "result": 8},
        {"num1": 9, "op": "-", "num2": 2, "result": 7},
    ]

    print("The hidden equations are:")
    for i, eq in enumerate(equations, 1):
        # In Python, ** is used for powers. For clarity, we'll show it as ^.
        # For division, we ensure integer output for clarity.
        if eq['op'] == '**':
            print(f"Line {i}: {eq['num1']} ^ {eq['num2']} = {eq['result']}")
        elif eq['op'] == '/':
            print(f"Line {i}: {eq['num1']} / {int(eq['num2'])} = {int(eq['result'])}")
        else:
            print(f"Line {i}: {eq['num1']} {eq['op']} {eq['num2']} = {eq['result']}")

solve_zoo_puzzle()