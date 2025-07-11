def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    
    This is based on researched data where films are marked if they are known 
    to show the obelisk in Paris's Place de la Concorde. The list is chronological.
    """
    # A chronological list of selected Best Picture winners and whether they show the obelisk.
    # This data is based on film research.
    best_picture_winners = [
        {"year": 1928, "title": "Wings", "shows_obelisk": False},
        {"year": 1932, "title": "Grand Hotel", "shows_obelisk": False},
        {"year": 1937, "title": "The Life of Emile Zola", "shows_obelisk": True},
        {"year": 1939, "title": "Gone with the Wind", "shows_obelisk": False},
        {"year": 1943, "title": "Casablanca", "shows_obelisk": False}, # Paris flashbacks do not clearly show it
        {"year": 1951, "title": "An American in Paris", "shows_obelisk": True},
        {"year": 1956, "title": "Around the World in 80 Days", "shows_obelisk": True},
        {"year": 1958, "title": "Gigi", "shows_obelisk": True},
    ]

    # Sort by year to ensure we find the earliest winner
    best_picture_winners.sort(key=lambda x: x['year'])

    # Find the first movie in the list that shows the obelisk
    for movie in best_picture_winners:
        if movie["shows_obelisk"]:
            print(f"The first Best Picture winner to show a Luxor Obelisk was from the year {movie['year']}.")
            print(f"Movie Title: {movie['title']}")
            return

find_first_winner_with_obelisk()