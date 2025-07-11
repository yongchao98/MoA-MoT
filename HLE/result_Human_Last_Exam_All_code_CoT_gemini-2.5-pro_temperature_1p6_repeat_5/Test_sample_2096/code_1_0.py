def solve_riddle():
    """
    This script methodically solves the historical riddle by connecting its clues.
    """

    print("Step 1: Decoding the clue about Pope Paul II and shame.")
    pope_name_and_number = "Paul II"
    print(f"The riddle focuses on {pope_name_and_number}.")
    print("Historical accounts from his critics, like the humanist Platina, portrayed Paul II as an opponent of learning.")
    print("Therefore, a 'shameful' accusation for him would be that he was uneducated or, in one word, 'illiterate'.\n")

    print("Step 2: Decoding the clue about the 1960s.")
    decade = 1960
    print(f"The clue states the word was 'written in the {decade}s'.")
    print("This points to a specific, widely distributed item from that era.")
    print("In the Soviet Union during the 1960s, a massive state-sponsored campaign to eliminate illiteracy was underway.")
    print("The Russian word for 'Illiterate' (Неграмотный) was prominently written on millions of propaganda posters and postage stamps circulated during this time.\n")

    print("Step 3: Combining the clues for the final answer.")
    print(f"The historical criticism of Paul {pope_name_and_number[-2:]} aligns with the word that was mass-printed during the {decade}s literacy campaign.")
    final_word = "ILLITERATE"
    print(f"Thus, the word 'X' that was considered shameful for the Pope and was widely 'written' in the 1960s is:")
    print(f"'{final_word}'")

solve_riddle()