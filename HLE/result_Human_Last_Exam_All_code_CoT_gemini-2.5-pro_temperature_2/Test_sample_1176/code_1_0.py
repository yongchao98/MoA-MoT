def solve_geology_task():
    """
    Solves the geology classification task by identifying each statement
    as consensus (C) or debated (D) and printing the results.
    """
    # Mapping statement number to its classification (C/D).
    # Based on an evaluation of the current state of academic consensus.
    classifications = {
        1: "D",
        2: "D",
        3: "D",
        4: "D",
        5: "C",
        6: "C",
        7: "C",
        8: "D",
        9: "D",
        10: "D",
    }

    print("Classification for each statement:")
    final_string = ""
    # Iterate through the sorted items to ensure order from 1 to 10
    for number, result in sorted(classifications.items()):
        print(f"Statement ({number}): {result}")
        final_string += result

    print("\n" + "="*20)
    print("Final combined string:")
    print(final_string)
    print("="*20)

solve_geology_task()