def solve_vienna_poet_riddle():
    """
    This function explains and solves the riddle about the English poet's
    surname found in a description of Vienna's boulevards.
    """

    # Step 1: The key geographical feature in Vienna is its famous
    # Ringstrasse boulevard, often just called "The Ring".
    location_clue = "The Ring (Ringstrasse)"
    
    # Step 2: In Russian, as per the puzzle's constraint, the word for "ring"
    # is "кольцо".
    russian_word = "кольцо"
    russian_pronunciation = "kol'tso"
    
    # Step 3: We need to find an English poet whose surname sounds like "кольцо".
    # The surname "Coleridge" is a very close match when transliterated into Russian.
    poet_surname_english = "Coleridge"
    poet_surname_russian = "Кольридж"
    poet_surname_russian_pronunciation = "Kol'ridzh"

    # Step 4: Present the logic and the final answer.
    print("Solving the riddle step-by-step:")
    print(f"1. The 'wide boulevards of Vienna' hints at the famous '{location_clue}'.")
    print(f"2. In the 'Russian translation', the word for 'ring' is '{russian_word}' (pronounced '{russian_pronunciation}').")
    print(f"3. An English poet's surname sounds very similar to this Russian word.")
    print(f"4. The poet is Samuel Taylor Coleridge. His name in Russian is '{poet_surname_russian}' (pronounced '{poet_surname_russian_pronunciation}').")
    print("\n" + "="*40)
    print("The English surname that appears in the description is:")
    print(f"'{poet_surname_english}'")
    print("="*40)

# Run the solver
solve_vienna_poet_riddle()