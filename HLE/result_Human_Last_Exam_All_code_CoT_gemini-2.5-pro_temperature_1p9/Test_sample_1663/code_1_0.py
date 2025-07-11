def solve_literary_puzzle():
    """
    This function provides the answer to the literary question about the Spanish poem.
    """
    poet = "Rafael Alberti"
    collection = "Sobre los ángeles (Concerning the Angels)"
    year = 1929
    poem_title = "Ángeles muertos (Dead Angels)"
    artwork_context = (
        "Rafael Alberti, who was also a painter, created a series of surrealist drawings "
        "to accompany the poems in his book. The verses are directed at the nightmarish, "
        "inescapable world depicted in both his poetry and his own illustrations for the book."
    )

    answer = f"""
The poet who wrote these verses is {poet}.
The lines are from the poem "{poem_title}", part of his collection "{collection}" published in {year}.

The work of art to which they were directed is not a separate piece by another artist; rather, {artwork_context}
"""
    print(answer)

solve_literary_puzzle()