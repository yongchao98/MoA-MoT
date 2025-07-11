def solve_riddle():
    """
    This script solves the riddle by laying out the logical steps connecting the clues.
    """
    
    # Clues from the riddle
    subject = "Pope Paul II"
    context = "Protecting European culture after the fall of Constantinople"
    shameful_quality = "Being uncultured or destructive of culture"
    
    # The solution word based on the historical context
    solution_word = "Barbarian"
    
    # The key to connecting the word to the final clue
    key_year = 1961
    key_work = "The poem 'Waiting for the Barbarians' by C. P. Cavafy"
    
    print("Step 1: Analyzing the historical context.")
    print(f"A Pope like {subject}, working to preserve culture, would consider it shameful to be a '{solution_word}'.")
    print("-" * 20)
    
    print("Step 2: Connecting the word to the 1960s clue.")
    print("The riddle states the word was 'written in the 1960s'.")
    print(f"This refers to {key_work}.")
    print(f"While the poem was composed earlier, a highly influential English translation was published in the year {key_year}.")
    print("This publication made the work famous in the English-speaking world during the 1960s.")
    print("-" * 20)
    
    print("Final Answer: The word X is:")
    print(solution_word)

solve_riddle()