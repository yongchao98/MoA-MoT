def find_first_best_picture_with_obelisk():
    """
    Identifies the first Best Picture winner to depict a Luxor Obelisk.

    This function simulates the research process by using a pre-compiled list of
    films relevant to Paris, Egypt, and the Best Picture award. It iterates
    chronologically to find the first movie that matches the criteria.
    """

    # Data on relevant films. The two Luxor Obelisks are in Luxor, Egypt, and Paris, France.
    # This list includes Best Picture winners set in Paris, as well as other notable films for context.
    # Each dictionary contains the title, award year, a boolean for winning Best Picture,
    # and a boolean for depicting a Luxor Obelisk.
    films_data = [
        {
            "title": "The Life of Emile Zola",
            "year": 1937,
            "won_best_picture": True,
            "depicts_obelisk": False  # Though set in Paris, it's not known for depicting the landmark.
        },
        {
            "title": "An American in Paris",
            "year": 1951,
            "won_best_picture": True,
            "depicts_obelisk": True  # Famously depicts Place de la Concorde in its ballet sequence.
        },
        {
            "title": "Around the World in 80 Days",
            "year": 1956,
            "won_best_picture": True,
            "depicts_obelisk": True  # Features a scene in Paris at the Place de la Concorde.
        },
        {
            "title": "Gigi",
            "year": 1958,
            "won_best_picture": True,
            "depicts_obelisk": True  # Features many Parisian landmarks.
        },
        {
            "title": "Cleopatra",
            "year": 1963,
            "won_best_picture": False, # Nominated, but did not win.
            "depicts_obelisk": True # Depicts ancient Egypt with obelisks.
        }
    ]

    # Sort films by year to ensure we find the earliest one first.
    films_data.sort(key=lambda x: x["year"])

    # Find the first film that meets both conditions.
    for film in films_data:
        if film["won_best_picture"] and film["depicts_obelisk"]:
            winner_title = film["title"]
            winner_year = film["year"]
            
            # The award year (e.g., 1951) refers to the year of film release.
            # The actual ceremony was in the following year (1952).
            # The 24th Academy Awards honored films of 1951.
            academy_award_number = winner_year - 1927
            
            print(f"The first Best Picture winner to depict a Luxor Obelisk was '{winner_title}'.")
            print(f"It won the {academy_award_number}th Academy Award for films released in the year {winner_year}.")
            
            # Stop after finding the first one.
            return

find_first_best_picture_with_obelisk()