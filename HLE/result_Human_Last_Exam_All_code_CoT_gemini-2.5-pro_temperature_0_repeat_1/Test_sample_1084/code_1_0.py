def solve():
    """
    This function checks the validity of six statements based on the derived condition q >= 1 + p.
    """
    statements = [
        # p, q
        (0.5, 2.0),
        (0.99, 2.0),
        (1.0, 2.0),
        (1.5, 2.0),
        (1.3, 2.6),
        (1.0, 1.5),
    ]

    result = []
    print("Checking the condition q >= 1 + p for each statement:")
    for i, (p, q) in enumerate(statements):
        is_true = q >= 1 + p
        char_result = 'Y' if is_true else 'N'
        result.append(char_result)
        print(f"Statement {i+1}: p={p}, q={q}. Condition: {q} >= 1 + {p} = {1+p}. This is {is_true}. -> {char_result}")

    final_answer = "".join(result)
    print("\nFinal answer string:")
    print(final_answer)
    
solve()