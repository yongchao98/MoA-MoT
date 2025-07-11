def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    The primary candidate location for this is the Place de la Concorde in Paris.
    """

    # A dictionary of early Best Picture winners, mapping the award ceremony year to the film title.
    # The data is ordered chronologically by the award ceremony.
    best_picture_winners = {
        1929: "The Broadway Melody",
        1930: "All Quiet on the Western Front",
        1931: "Cimarron",
        1932: "Grand Hotel",
        1934: "Cavalcade",
        1935: "It Happened One Night",
        1936: "Mutiny on the Bounty",
        1937: "The Great Ziegfeld",
        1938: "The Life of Emile Zola",
        1939: "You Can't Take It with You",
        1940: "Gone with the Wind",
        1941: "Rebecca",
        1942: "How Green Was My Valley",
        1943: "Mrs. Miniver",
        1944: "Casablanca",
        1945: "Going My Way",
        1946: "The Lost Weekend",
        1947: "The Best Years of Our Lives",
        1948: "Gentleman's Agreement",
        1949: "Hamlet",
        1950: "All the King's Men",
        1951: "All About Eve",
        1952: "An American in Paris",  # Set in Paris
        1953: "The Greatest Show on Earth",
        1954: "From Here to Eternity",
        1955: "On the Waterfront",
        1956: "Marty",
        1957: "Around the World in 80 Days", # Also set in Paris
    }

    # A list of films known to prominently feature the Luxor Obelisk in Paris.
    films_with_paris_obelisk = [
        "An American in Paris",
        "Around the World in 80 Days"
    ]

    # Iterate through the winners chronologically by sorting the award years.
    for year in sorted(best_picture_winners.keys()):
        winner_title = best_picture_winners[year]
        
        # Check if the winner is one of the films known to feature the obelisk.
        if winner_title in films_with_paris_obelisk:
            print(f"The first Best Picture winner to depict a Luxor Obelisk on-screen is: {winner_title}")
            return

find_first_winner_with_obelisk()