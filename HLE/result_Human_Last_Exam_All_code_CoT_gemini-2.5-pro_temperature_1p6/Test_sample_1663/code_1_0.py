import textwrap

def present_answer():
    """
    This function stores and prints the information about the poem and artwork.
    """
    poet = "Rafael Alberti"
    artwork = "Guernica"
    artist = "Pablo Picasso"
    poem_collection = "A la pintura (To Painting)"

    verses = [
        "No one leaves from here. Nobody.",
        "Neither the mystic nor the suicidal.",
        "And it's useless,",
        "All escape is useless",
        "(Not even from below",
        " or from above)."
    ]

    explanation = (
        f"These verses were written by the Spanish poet {poet}. They were directed "
        f"at the masterpiece painting '{artwork}' by {artist}. The poem is "
        f"part of Alberti's larger work, '{poem_collection}', which is a tribute "
        f"to the art of painting."
    )

    print("--- Poetry and Art Inquiry ---")
    print("\n[The Verses]")
    for line in verses:
        print(f'  "{line}"')
    
    print("\n[Analysis]")
    # The textwrap module helps format the paragraph nicely in the terminal.
    for line in textwrap.wrap(explanation, width=70):
        print(line)

# Execute the function to display the answer.
present_answer()
