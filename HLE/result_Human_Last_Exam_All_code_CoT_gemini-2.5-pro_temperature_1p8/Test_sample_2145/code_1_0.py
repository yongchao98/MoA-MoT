def solve_trivia_puzzle():
    """
    Solves a multi-part trivia puzzle to find a hidden word.
    """

    # --- Question 1 ---
    answer1 = "Slumdog Millionaire"
    letter1 = answer1[0]
    explanation1 = "The sequence A, A, D, A, C, A, D, D, A represents the answers given by the protagonist in the Oscar-winning film 'Slumdog Millionaire'."

    # --- Question 2 ---
    answer2 = "Enemies"
    letter2 = answer2[0]
    explanation2 = "Antagonists or 'Enemies' were placed on the left in Soviet cinema and are depicted at the bottom of historic 'Last Judgment' diagrams."

    # --- Question 3 ---
    answer3 = "Alf"
    letter3 = answer3[0]
    explanation3 = "The jokes about eating cats ('Cheshire Salad', 'Siamese Hotpot') refer to the title character of the TV series 'Alf'."

    # --- Question 4 ---
    answer4 = "Maugli"
    letter4 = answer4[0]
    explanation4 = "The red dogs in the Soviet-era cartoon 'Maugli' were a caricature of Mao's Red Guards; Mao came to power through revolution."

    # --- Combine and Print ---
    hidden_word = letter1 + letter2 + letter3 + letter4

    print("Here are the answers to the trivia questions and the hidden word:")
    print("-" * 60)
    print(f"(1) {explanation1}\n    The answer is '{answer1}', so the first letter is '{letter1}'.\n")
    print(f"(2) {explanation2}\n    The answer is '{answer2}', so the first letter is '{letter2}'.\n")
    print(f"(3) {explanation3}\n    The answer is '{answer3}', so the first letter is '{letter3}'.\n")
    print(f"(4) {explanation4}\n    The answer is '{answer4}', so the first letter is '{letter4}'.\n")
    print("-" * 60)
    print("Combining the first letters gives the final 'equation':")
    print(f"'{letter1}' + '{letter2}' + '{letter3}' + '{letter4}' = {hidden_word}")

solve_trivia_puzzle()

<<<SEAM>>>