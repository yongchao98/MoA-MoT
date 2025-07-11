def solve_riddle():
    """
    This script solves the riddle by breaking it down into two parts
    and connecting them with a pun.
    """

    # Part 1: Analyze the clue about mercury salts.
    # The clue "Another X Y has almost ceased to be used due to the toxicity of mercury salts"
    # refers to the historical process of making felt. Mercuric nitrate was used to make
    # felt for hats, a practice now discontinued due to its toxicity.
    # This identifies "X Y" as "Felt Hat".
    X = "Felt"
    Y = "Hat"

    # Part 2: Connect "Felt Hat" to the image.
    # The image shows hockey player Nikita Kucherov at the 2019 NHL Awards.
    # At this event, he won the "Hart" Memorial Trophy.
    # The word "Hart" sounds the same as "Hat".
    connection = "The 'Hart' Trophy sounds like 'Hat'."

    # Conclusion: The two clues connect through this pun.
    final_answer = f"{X} {Y}"

    print(f"Clue 1 Analysis: The image shows the winner of the 'Hart' Memorial Trophy.")
    print(f"Clue 2 Analysis: The process of making '{X}' for hats with mercury salts has ceased.")
    print(f"The Connection (Pun): 'Hart' sounds like '{Y}'.")
    print("-" * 20)
    print(f"Therefore, X Y is: {final_answer}")


solve_riddle()
<<<Felt Hat>>>