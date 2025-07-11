def find_first_winner_with_obelisk():
    """
    Finds the first Best Picture winner to depict a Luxor Obelisk by
    cross-referencing a list of winners with a list of films known
    to feature the landmark.
    """

    # Data Source 1: A dictionary of Academy Award for Best Picture winners.
    # The key is the year of the film's release.
    best_picture_winners = {
        1928: "Wings",
        1929: "The Broadway Melody",
        1930: "All Quiet on the Western Front",
        1931: "Cimarron",
        1932: "Grand Hotel",
        1933: "Cavalcade",
        1934: "It Happened One Night",
        1935: "Mutiny on the Bounty",
        1936: "The Great Ziegfeld",
        1937: "The Life of Emile Zola",
        1938: "You Can't Take It with You",
        1939: "Gone with the Wind",
        1940: "Rebecca",
        1941: "How Green Was My Valley",
        1942: "Mrs. Miniver",
        1943: "Casablanca",
        1944: "Going My Way",
        1945: "The Lost Weekend",
        1946: "The Best Years of Our Lives",
        1947: "Gentleman's Agreement",
        1948: "Hamlet",
        1949: "All the King's Men",
        1950: "All About Eve",
        1951: "An American in Paris",
        1952: "The Greatest Show on Earth",
        1953: "From Here to Eternity",
        1954: "On the Waterfront",
        1955: "Marty",
        1956: "Around the World in 80 Days",
        1957: "The Bridge on the River Kwai",
        1958: "Gigi",
        1959: "Ben-Hur"
    }

    # Data Source 2: A list of films known to depict one of the two Luxor Obelisks.
    films_with_luxor_obelisk = [
        "An American in Paris",        # Features the obelisk in Place de la Concorde, Paris
        "Around the World in 80 Days", # Features the obelisk in Place de la Concorde, Paris
        "Gigi",                        # Features the obelisk in Place de la Concorde, Paris
        "Death on the Nile"            # Features the obelisk at Luxor Temple, Egypt
    ]

    # Sort the dictionary keys (years) to process the films chronologically.
    sorted_years = sorted(best_picture_winners.keys())

    # Iterate through the years and find the first winner that is in our obelisk list.
    for year in sorted_years:
        title = best_picture_winners[year]
        if title in films_with_luxor_obelisk:
            print(f"The first Academy Award winner for Best Picture to depict a Luxor Obelisk is the winner from {year}:")
            print(title)
            return

# Run the function to find and print the answer.
find_first_winner_with_obelisk()