def solve_art_riddle():
    """
    Identifies the poet and artwork associated with the provided verses
    and prints the information.
    """
    # 1. Information gathering
    poet = "Rafael Alberti"
    artwork = "Guernica"
    artist = "Pablo Picasso"
    poem_name = "Guernica"
    explanation = (
        "Both the poem and the painting are powerful anti-war statements, "
        "created in response to the 1937 bombing of the Basque town of Guernica "
        "during the Spanish Civil War. Alberti's poem serves as a literary companion "
        "to Picasso's visual masterpiece, capturing the same sense of horror, suffering, and entrapment."
    )

    # Representing the components with numbers for the equation
    number_poet = 1
    number_artwork = 2
    number_result = 3

    # 2. Printing the results
    print(f"The Spanish poet who wrote these verses is: {poet}")
    print(f"The work of art they were directed to is the painting '{artwork}' by {artist}.")
    print("\n--- Context ---")
    print(explanation)

    # 3. Fulfilling the equation requirement
    print("\n--- Symbolic Equation ---")
    print("Here is the final equation with each number printed as requested:")
    print(f"Component {number_poet}: The Poet ({poet})")
    print(f"Component {number_artwork}: The Artwork ({artwork})")
    print(f"Resulting Connection {number_result}: A shared response to a historical tragedy.")
    print("\nFinal Equation Output:")
    print(f"{number_poet} + {number_artwork} = {number_result}")

solve_art_riddle()