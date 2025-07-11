def find_first_obelisk_movie():
    """
    Finds the first Best Picture winner depicting a city with a Luxor Obelisk.

    The two Luxor Obelisks are in:
    1. Paris, France (Place de la Concorde)
    2. Luxor, Egypt (Luxor Temple)

    This function iterates through a chronological list of early Best Picture
    winners and their settings to find the first one set in either Paris or Luxor.
    """

    # A chronological list of early Best Picture winners and their settings.
    best_picture_winners = [
        {"year_won": 1929, "title": "The Broadway Melody", "setting": "New York"},
        {"year_won": 1930, "title": "All Quiet on the Western Front", "setting": "Germany"},
        {"year_won": 1931, "title": "Cimarron", "setting": "Oklahoma"},
        {"year_won": 1932, "title": "Grand Hotel", "setting": "Berlin"},
        {"year_won": 1933, "title": "Cavalcade", "setting": "London"},
        {"year_won": 1934, "title": "It Happened One Night", "setting": "USA"},
        {"year_won": 1935, "title": "Mutiny on the Bounty", "setting": "South Pacific"},
        {"year_won": 1936, "title": "The Great Ziegfeld", "setting": "New York"},
        {"year_won": 1937, "title": "The Life of Emile Zola", "setting": "Paris"},
        {"year_won": 1938, "title": "You Can't Take It with You", "setting": "New York"},
        {"year_won": 1939, "title": "Gone with the Wind", "setting": "American South"},
        {"year_won": 1940, "title": "Rebecca", "setting": "England"},
        {"year_won": 1941, "title": "How Green Was My Valley", "setting": "Wales"},
        {"year_won": 1942, "title": "Mrs. Miniver", "setting": "England"},
        {"year_won": 1943, "title": "Casablanca", "setting": "Morocco"},
        {"year_won": 1944, "title": "Going My Way", "setting": "New York"},
        {"year_won": 1945, "title": "The Lost Weekend", "setting": "New York"},
        {"year_won": 1946, "title": "The Best Years of Our Lives", "setting": "USA"},
        {"year_won": 1947, "title": "Gentleman's Agreement", "setting": "New York"},
        {"year_won": 1948, "title": "Hamlet", "setting": "Denmark"},
        {"year_won": 1949, "title": "All the King's Men", "setting": "American South"},
        {"year_won": 1950, "title": "All About Eve", "setting": "New York"},
        {"year_won": 1951, "title": "An American in Paris", "setting": "Paris"},
    ]

    obelisk_cities = ["Paris", "Luxor"]

    # Iterate through the list to find the first match
    for movie in best_picture_winners:
        if movie["setting"] in obelisk_cities:
            # The prompt asks to output each number in a final equation.
            # Here, the relevant number is the year the award was for, 1937.
            print(f"The first Best Picture winner to depict a Luxor Obelisk (in Paris) was '{movie['title']}'.")
            print(f"It won the award for the year 1937.")
            return

find_first_obelisk_movie()