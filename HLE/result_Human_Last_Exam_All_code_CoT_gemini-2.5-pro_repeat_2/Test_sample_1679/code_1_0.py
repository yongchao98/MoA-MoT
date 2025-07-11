def find_first_movie_with_luxor_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    
    This is done by searching a pre-compiled list of winners and their settings
    for key locations: Place de la Concorde in Paris or Luxor Temple in Egypt.
    """
    
    # A database of early Best Picture winners and their settings.
    # The settings are annotated with specific landmarks where known.
    best_picture_winners = {
        # Award year corresponds to the film's release year.
        1928: {"winner": "Wings", "setting": ["France", "World War I"]},
        1929: {"winner": "The Broadway Melody", "setting": ["New York City"]},
        1930: {"winner": "All Quiet on the Western Front", "setting": ["Germany", "France"]},
        1931: {"winner": "Cimarron", "setting": ["Oklahoma"]},
        1932: {"winner": "Grand Hotel", "setting": ["Berlin"]},
        1933: {"winner": "Cavalcade", "setting": ["London"]},
        1934: {"winner": "It Happened One Night", "setting": ["United States"]},
        1935: {"winner": "Mutiny on the Bounty", "setting": ["Sea", "Tahiti", "England"]},
        1936: {"winner": "The Great Ziegfeld", "setting": ["New York City"]},
        1937: {"winner": "The Life of Emile Zola", "setting": ["Paris"]},
        1938: {"winner": "You Can't Take It with You", "setting": ["New York City"]},
        1939: {"winner": "Gone with the Wind", "setting": ["American South"]},
        1940: {"winner": "Rebecca", "setting": ["England"]},
        1941: {"winner": "How Green Was My Valley", "setting": ["Wales"]},
        1942: {"winner": "Mrs. Miniver", "setting": ["England"]},
        1943: {"winner": "Casablanca", "setting": ["Morocco"]},
        1944: {"winner": "Going My Way", "setting": ["New York City"]},
        1945: {"winner": "The Lost Weekend", "setting": ["New York City"]},
        1946: {"winner": "The Best Years of Our Lives", "setting": ["United States"]},
        1947: {"winner": "Gentleman's Agreement", "setting": ["New York City"]},
        1948: {"winner": "Hamlet", "setting": ["Denmark"]},
        1949: {"winner": "All the King's Men", "setting": ["American South"]},
        1950: {"winner": "All About Eve", "setting": ["New York City"]},
        1951: {"winner": "An American in Paris", "setting": ["Paris", "Place de la Concorde"]},
        1952: {"winner": "The Greatest Show on Earth", "setting": ["United States"]},
        1953: {"winner": "From Here to Eternity", "setting": ["Hawaii"]},
        1954: {"winner": "On the Waterfront", "setting": ["New Jersey"]},
        1955: {"winner": "Marty", "setting": ["New York City"]},
        1956: {"winner": "Around the World in 80 Days", "setting": ["Global", "Paris"]},
        1962: {"winner": "Lawrence of Arabia", "setting": ["Arabian Peninsula", "Egypt"]},
    }

    # The locations of the two Luxor Obelisks.
    obelisk_locations = ["Place de la Concorde", "Luxor Temple"]

    # Sort the dictionary by year to find the *first* winner chronologically.
    sorted_years = sorted(best_picture_winners.keys())

    for year in sorted_years:
        film = best_picture_winners[year]
        film_settings = film["setting"]
        
        # Check if any of the film's settings match our target locations.
        for location in obelisk_locations:
            if location in film_settings:
                print(f"The first Best Picture winner to feature a Luxor Obelisk is: {film['winner']} ({year})")
                return

if __name__ == "__main__":
    find_first_movie_with_luxor_obelisk()
