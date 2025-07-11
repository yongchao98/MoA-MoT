def solve_riddle():
    """
    This function explains the steps to solve the Vienna boulevards riddle.
    """
    print("The riddle describes the 'wide boulevards' of Vienna.")
    print("Step 1: Translate the key descriptive word 'wide' into Russian.")
    
    russian_word = "широкий"
    transliteration = "shirokiy"
    print(f"The Russian word for 'wide' is '{russian_word}', which is pronounced roughly as '{transliteration}'.")
    
    print("\nStep 2: Find the English poet's surname that sounds like this Russian word.")
    
    poet_surname = "Shelley"
    print(f"The pronunciation '{transliteration}' is phonetically very similar to the surname of the English Romantic poet, Percy Bysshe {poet_surname}.")
    
    print("\n--- The Final 'Equation' ---")
    print(f"Russian for 'wide' ('{russian_word}') ==> sounds like ==> English poet's name ('{poet_surname}')")

solve_riddle()