def find_movie_with_obelisk():
    """
    This function identifies the first Best Picture Oscar winner to depict a Luxor Obelisk.
    """

    # There are two Luxor Obelisks.
    obelisk_location_1 = "Luxor Temple, Luxor, Egypt"
    obelisk_location_2 = "Place de la Concorde, Paris, France"

    # After reviewing the list of Best Picture winners, the first one
    # to be set in a location with a Luxor Obelisk is "An American in Paris".
    movie_title = "An American in Paris"
    release_year = 1951
    award_ceremony_number = 24
    award_year = 1952 # The ceremony for 1951 films was in 1952.

    # The film features the obelisk in the Place de la Concorde during its famous ballet sequence.
    depicted_obelisk_location = obelisk_location_2

    print("Finding the first Best Picture winner to show a Luxor Obelisk:")
    print("-" * 60)
    print(f"The two Luxor Obelisks are located in:\n1. {obelisk_location_1}\n2. {obelisk_location_2}")
    print("-" * 60)
    print(f"The first Academy Award winner for Best Picture to depict one of these is:")
    print(f"Movie Title: {movie_title}")
    print(f"Release Year: {release_year}")
    print(f"It won the {award_ceremony_number}rd Academy Award for Best Picture.")
    print(f"The obelisk shown is the one in the {depicted_obelisk_location}.")
    print("-" * 60)

find_movie_with_obelisk()