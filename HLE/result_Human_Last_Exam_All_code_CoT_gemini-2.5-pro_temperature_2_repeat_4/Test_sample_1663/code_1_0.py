def find_poet_and_art():
    """
    This function identifies the poet and the artwork associated with the verses.
    """
    poet = "Rafael Alberti"
    poem_info = "The poem is 'Ángeles muertos' ('Dead Angels') from his 1929 collection 'Sobre los ángeles' ('Concerning the Angels')."
    
    artist = "El Greco"
    artwork = "The Burial of the Count of Orgaz (El entierro del conde de Orgaz)"

    print("The poet who wrote these verses is:")
    print(f"- {poet}")
    print("\n" + poem_info)
    print("\nThe poem was inspired by and directed to a work of art by:")
    print(f"- Artist: {artist}")
    print(f"- Artwork: {artwork}")

if __name__ == "__main__":
    find_poet_and_art()