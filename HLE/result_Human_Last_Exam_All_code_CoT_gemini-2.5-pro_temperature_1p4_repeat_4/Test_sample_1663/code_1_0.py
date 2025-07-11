def find_poet_and_artwork():
    """
    This function stores and prints the information about the poet
    and the artwork associated with the provided verses.
    """
    poet = "Rafael Alberti"
    artwork_subject = "The 'Black Paintings' (Pinturas Negras)"
    artist = "Francisco de Goya"
    poem_collection = "A la pintura (To Painting)"
    
    # The verses are from a poem dedicated to the color black,
    # inspired by Goya's dark and intense late works.

    print("The Spanish poet who wrote these verses is:")
    print(f"  - {poet}")
    print("\nThe verses were directed to a famous series of paintings:")
    print(f"  - Artist: {artist}")
    print(f"  - Work of Art: {artwork_subject}")
    print(f"\nThe poem is from his collection '{poem_collection}'.")

if __name__ == "__main__":
    find_poet_and_artwork()