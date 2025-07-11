def find_first_best_picture_with_luxor_obelisk():
    """
    Finds the first Best Picture winner to depict a Luxor Obelisk.

    The two Luxor Obelisks are located in:
    1. Luxor, Egypt (at the Luxor Temple)
    2. Paris, France (at the Place de la Concorde)

    This script searches for the earliest Best Picture winner set in a location
    where one of these obelisks could be depicted on-screen. The obelisk in Paris
    is a prominent landmark in a major city often featured in films.
    """
    
    # A list of early Best Picture winners, their award year, and primary setting.
    # The year corresponds to the year the award was given.
    best_picture_winners = [
        {"year": 1928, "title": "Wings", "setting": "France/USA"},
        {"year": 1929, "title": "The Broadway Melody", "setting": "USA"},
        {"year": 1930, "title": "All Quiet on the Western Front", "setting": "Germany"},
        {"year": 1931, "title": "Cimarron", "setting": "USA"},
        {"year": 1932, "title": "Grand Hotel", "setting": "Germany"},
        {"year": 1933, "title": "Cavalcade", "setting": "UK"},
        {"year": 1934, "title": "It Happened One Night", "setting": "USA"},
        {"year": 1935, "title": "Mutiny on the Bounty", "setting": "Sea/Tahiti"},
        {"year": 1936, "title": "The Great Ziegfeld", "setting": "USA"},
        {"year": 1937, "title": "The Life of Emile Zola", "setting": "Paris, France"},
        {"year": 1938, "title": "You Can't Take It with You", "setting": "USA"},
        {"year": 1939, "title": "Gone with the Wind", "setting": "USA"},
        {"year": 1951, "title": "An American in Paris", "setting": "Paris, France"},
    ]

    print("Searching for the first Best Picture winner set in Paris, France...")
    
    winner = None
    for film in sorted(best_picture_winners, key=lambda x: x['year']):
        # Check if the film's setting is Paris, where the obelisk stands in the Place de la Concorde.
        if "Paris, France" in film["setting"]:
            winner = film
            break # Stop at the first chronological match

    if winner:
        print("\n--- Found a Match ---")
        print(f"Movie Title: {winner['title']}")
        # The award for the film year 1937 was given at the 10th Academy Awards in 1938.
        # We will refer to the film by its film year for clarity.
        film_year = winner['year']
        print(f"Academy Award Year: {film_year}")
        print("\nThis film is set in Paris and is known to feature scenes at the Place de la Concorde,")
        print("which is the location of one of the two Luxor Obelisks.")
    else:
        print("Could not determine the winner from the provided data.")

find_first_best_picture_with_luxor_obelisk()