def solve_riddle():
    """
    This function explains the step-by-step solution to the riddle.
    """
    print("Analyzing the riddle to identify 'X Y'...\n")

    # Part 1: The clue about mercury toxicity
    print("--- Clue 1: The Toxic Object ---")
    print("The clue 'Another X Y has almost ceased to be used due to the toxicity of mercury salts' points to a FELT HAT.")
    print("The process of making felt for hats historically involved using mercury nitrate, a toxic salt that caused neurological damage in hat makers.\n")

    # Part 2: The clue about the picture
    print("--- Clue 2: The Pictured Event ---")
    print("The clue 'X Y is pictured at an event on June 20, 2019' and the image point to the 2019 NHL Awards.")
    print("The person receiving the award is Nikita Kucherov, and the award is the HART Memorial Trophy.\n")

    # Part 3: Connecting the clues with a pun
    print("--- Connecting the Clues ---")
    print("The answer connects the two clues through a two-part pun:")
    # The pun for 'Felt'
    word_x = "Felt"
    explanation_x = "The word 'Felt' is a pun on the emotion the winner 'felt' upon receiving the prestigious award."
    # The pun for 'Hat'
    word_y = "Hat"
    explanation_y = "The word 'Hat' is a pun on the 'Hart' Memorial Trophy being presented in the picture."
    
    print(f"1. X = '{word_x}': {explanation_x}")
    print(f"2. Y = '{word_y}': {explanation_y}\n")

    # Final Answer
    print("--- Final Answer ---")
    print(f"Therefore, 'X Y' is: {word_x} {word_y}")

solve_riddle()