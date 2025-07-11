def find_poet_and_artwork():
    """
    This function identifies and prints the author and subject of the provided verses.
    """
    poet = "Rafael Alberti"
    artwork = "Pablo Picasso's painting, 'Guernica'"
    poem_collection = "'A la pintura' (To Painting)"

    print("The provided verses were written by the Spanish poet:")
    print(f"- {poet}")
    print("\nThey were directed to the work of art:")
    print(f"- {artwork}")
    print(f"\nThe verses are part of a section dedicated to the painting within Alberti's larger collection, {poem_collection}.")

# Execute the function to display the information
find_poet_and_artwork()