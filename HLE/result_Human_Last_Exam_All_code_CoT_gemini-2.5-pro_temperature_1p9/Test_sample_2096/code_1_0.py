def solve_riddle():
    """
    This script solves the lateral thinking puzzle step-by-step.
    """
    print("Thinking Process to solve the riddle:")
    print("-" * 40)

    # Step 1: Analyze the first clue about the Pope.
    print("Step 1: Analyze the clue 'being X was considered shameful for the Pope.'")
    pope_clue_analysis = (
        "Pope Paul II was a Renaissance figure. During this 'rebirth' of learning,\n"
        "for a Pope to be unable to read or write would be a great source of shame.\n"
        "This suggests the word X could be 'ILLITERATE'."
    )
    print(pope_clue_analysis)
    print("")

    # Step 2: Analyze the second clue about the 1960s.
    print("Step 2: Analyze the clue '\"X\" was written in the 1960s.'")
    sixties_clue_analysis = (
        "This is a play on words. It doesn't mean the word was invented in the 1960s.\n"
        "It refers to the word literally being written down. In the 1960s, a famous piece of graffiti\n"
        "in London stated 'Clapton is God'. A witty response was written underneath: '...and he's illiterate'.\n"
        "This act of writing fits the clue perfectly."
    )
    print(sixties_clue_analysis)
    print("")

    # Step 3: Combine the clues to find the answer.
    print("Step 3: Combine the analysis of both clues.")
    final_answer = "ILLITERATE"
    conclusion = (
        f"The word that satisfies both conditions is '{final_answer}'.\n"
        f"- Being {final_answer} was shameful for a Pope.\n"
        f"- The word '{final_answer}' was famously written as graffiti in the 1960s."
    )
    print(conclusion)
    print("-" * 40)

    print(f"\nThe final answer is the word: {final_answer}")

solve_riddle()