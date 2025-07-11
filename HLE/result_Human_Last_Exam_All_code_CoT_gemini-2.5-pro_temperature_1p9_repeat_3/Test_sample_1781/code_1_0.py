def solve_history_puzzle():
    """
    This function analyzes the historical question and determines the best answer
    from the given choices.
    """
    # The question concerns King Philip II of France and has two parts:
    # 1. The year the royal title shifted from "King of the Franks" (personal) to "King of France" (territorial).
    #    This event began in 1190.
    # 2. The author of the biography that explains the king's epithet, "Augustus".
    #    This author was the monk Rigord.

    # The ideal answer is (1190, Rigord). Now let's evaluate the options provided.
    options = {
        'A': (1190, 'Suetonius'),
        'B': (1250, 'Joinville'),
        'C': (1223, 'Rigord'),
        'D': (1789, 'Voltaire'),
        'E': (1190, 'Baldwin')
    }

    # Analysis of Option C:
    # Author: Rigord is the correct biographer for Philip II and the source of his epithet "Augustus".
    # Year: 1223 is the year of Philip II's death. This date represents the culmination of his reign,
    # during which the "morphing" of the monarchy was completed and solidified.
    # While 1190 is the start date, pairing the correct author with the end date of the reign
    # makes this the strongest and most historically coherent choice among the options.

    correct_choice = 'C'
    year, author = options[correct_choice]

    print("Step 1: The question refers to the French monarch Philip II 'Augustus'.")
    print("Step 2: The shift in his title to reflect territoriality ('King of France') began in 1190.")
    print("Step 3: The source of his epithet 'Augustus' was his biographer, Rigord.")
    print("Step 4: Evaluating the choices, option C correctly identifies the biographer, Rigord.")
    print("Step 5: The year in option C is 1223, which is the year Philip II died, marking the end of the era containing the 'morphing'.")
    print("\nFinal Answer's components:")
    print(f"Year: {year}")
    print(f"Author: {author}")

solve_history_puzzle()