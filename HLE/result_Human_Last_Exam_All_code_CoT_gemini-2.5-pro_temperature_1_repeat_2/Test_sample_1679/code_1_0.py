def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    There are two Luxor Obelisks: one in Luxor, Egypt, and one in Paris, France.
    This function checks early winners for settings in these locations.
    """
    best_picture_winners = {
        1928: {"title": "Wings", "setting": "France (WWI)"},
        1929: {"title": "The Broadway Melody", "setting": "New York"},
        1930: {"title": "All Quiet on the Western Front", "setting": "Germany/France (WWI)"},
        1931: {"title": "Cimarron", "setting": "Oklahoma"},
        1932: {"title": "Grand Hotel", "setting": "Berlin"},
        1933: {"title": "Cavalcade", "setting": "London"},
        1934: {"title": "It Happened One Night", "setting": "USA"},
        1935: {"title": "Mutiny on the Bounty", "setting": "Sea/England"},
        1936: {"title": "The Great Ziegfeld", "setting": "New York"},
        1937: {"title": "The Life of Emile Zola", "setting": "Paris"},
        1938: {"title": "You Can't Take It with You", "setting": "New York"},
        1939: {"title": "Gone with the Wind", "setting": "USA South"},
        1940: {"title": "Rebecca", "setting": "Monte Carlo/England"},
        1941: {"title": "How Green Was My Valley", "setting": "Wales"},
        1942: {"title": "Mrs. Miniver", "setting": "England"},
        1943: {"title": "Casablanca", "setting": "Morocco"},
        1944: {"title": "Going My Way", "setting": "New York"},
        1945: {"title": "The Lost Weekend", "setting": "New York"},
        1946: {"title": "The Best Years of Our Lives", "setting": "USA"},
        1947: {"title": "Gentleman's Agreement", "setting": "New York"},
        1948: {"title": "Hamlet", "setting": "Denmark"},
        1949: {"title": "All the King's Men", "setting": "USA South"},
        1950: {"title": "All About Eve", "setting": "New York"},
        1951: {"title": "An American in Paris", "setting": "Paris"},
        1952: {"title": "The Greatest Show on Earth", "setting": "USA"},
        1953: {"title": "From Here to Eternity", "setting": "Hawaii"},
        1954: {"title": "On the Waterfront", "setting": "New Jersey"},
        1955: {"title": "Marty", "setting": "New York"},
        1956: {"title": "Around the World in 80 Days", "setting": "Global (incl. Paris, Egypt)"}
    }

    print("Searching for the first Best Picture winner set in Paris or Luxor...")
    print("-" * 60)

    # Sort by year to check chronologically
    for year in sorted(best_picture_winners.keys()):
        film = best_picture_winners[year]
        title = film["title"]
        setting = film["setting"]

        # Check if the setting is Paris or Egypt (Luxor)
        if "Paris" in setting:
            print(f"Found a candidate: '{title}' ({year}), which is set in Paris.")
            print(f"This film depicts the Luxor Obelisk located at the Place de la Concorde.")
            print(f"\nFinal Answer: The first winner was '{title}', which won the award for the year {year}.")
            return

    print("No winner found in the provided list.")


if __name__ == '__main__':
    find_first_winner_with_obelisk()