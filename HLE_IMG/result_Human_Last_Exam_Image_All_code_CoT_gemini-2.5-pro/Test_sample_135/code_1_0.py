def find_artist():
    """
    This function identifies the artist of the provided artwork detail.
    The identification is based on art historical analysis of the style,
    subject matter, and inscription.
    """
    # Information gathered through research
    artist = "Luca della Robbia"
    artwork_title = "Tabernacle of the Sacrament"
    detail_scene = "Institution of the Eucharist"
    location = "Church of Santa Maria a Peretola, Florence"
    date = "c. 1441-1443"

    # Print the findings
    print(f"The image shows a detail from the '{artwork_title}'.")
    print(f"The specific scene is the '{detail_scene}'.")
    print(f"This work is located in the {location} and dates to {date}.")
    print("\nThe author of this detail of the Altarpiece is:")
    print(artist)

# Execute the function to find and print the author's name.
find_artist()