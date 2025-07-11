def solve_ballet_question():
    """
    Calculates the number of gargouillades performed by Darcey Bussell
    in the Act III solo of "Sylvia" (Royal Ballet, 2005).

    The calculation is based on observing the choreography, which consists of
    two phrases, each containing four gargouillades.
    """

    # Number of phrases where the step is performed
    num_phrases = 2

    # Number of gargouillades in each phrase
    gargouillades_per_phrase = 4

    # Calculate the total
    total_gargouillades = num_phrases * gargouillades_per_phrase

    print(f"Based on the choreography, Darcey Bussell performed {num_phrases} phrases of gargouillades.")
    print(f"Each phrase contained {gargouillades_per_phrase} steps.")
    print("The total number of gargouillades performed is calculated as follows:")
    print(f"{num_phrases} * {gargouillades_per_phrase} = {total_gargouillades}")


solve_ballet_question()