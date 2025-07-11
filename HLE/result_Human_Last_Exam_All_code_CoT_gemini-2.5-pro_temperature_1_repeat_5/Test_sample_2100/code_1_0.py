def solve_riddle():
    """
    This function solves the historical riddle by identifying the person based on the clues.
    """
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    grave_clue = "Grave location unknown until the late 1980s"

    # Andrei Voznesensky wrote the story for "Juno and Avos," a famous Russian rock opera.
    # The opera is about the life and death of the Russian explorer Nikolai Rezanov.
    # Rezanov's grave was lost after being destroyed and was rediscovered/re-memorialized in the late 20th century.
    last_name = "Rezanov"

    print(f"The poet, {poet}, famously wrote about the Russian explorer Nikolai Rezanov.")
    print(f"This explorer fits the description of a 'sailor' and his grave's location was a mystery for a long time.")
    print(f"The last name of the man is: {last_name}")

solve_riddle()