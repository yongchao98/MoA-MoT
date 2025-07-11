def solve_bulgakov_puzzle():
    """
    This function identifies and prints the solution to the literary puzzle
    based on parallels in Bulgakov's "The Master and Margarita".
    """
    # Step 1: Identify the Moscow character tormented by a sparrow in Chapter 18.
    # In Chapter 18, "Unlucky Visitors," Andrei Fokich Sokov, the buffet-master,
    # is tormented by an "impudent and bold-faced sparrow" after leaving Woland's apartment.
    moscow_character = "Sokov"

    # Step 2 & 3: Identify the corresponding bird and character in Jerusalem.
    # The parallel scene is in Chapter 2, "Pontius Pilate." During the trial,
    # while Pilate is suffering from a terrible migraine, a swallow flies into the colonnade.
    # This bird is a barn swallow.
    jerusalem_bird = "barn swallow"
    
    # Step 4: Identify who the bird flies around.
    # The barn swallow flies around Pontius Pilate during his moment of great stress and torment.
    jerusalem_character = "Pontius Pilate"

    # Step 5: Format and print the final answer.
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"
    print(final_answer)

solve_bulgakov_puzzle()