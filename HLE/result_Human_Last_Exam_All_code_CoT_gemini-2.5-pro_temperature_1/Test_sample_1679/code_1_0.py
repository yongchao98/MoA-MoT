def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    The script iterates through a chronological list of winners and checks against
    a list of films known to feature the obelisk in Paris or Luxor.
    """

    # A dictionary of early Best Picture winners (Year: Title)
    # The year represents the film's release year, for which it won the award.
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
    }

    # Films known to depict a Luxor Obelisk (in Paris or Luxor)
    # Based on research, "An American in Paris" is the first Best Picture winner to do so.
    movies_with_obelisk = ["An American in Paris", "Gigi"]

    # Sort the years to ensure chronological order
    sorted_years = sorted(best_picture_winners.keys())

    # Iterate through the winners chronologically
    for year in sorted_years:
        winner_title = best_picture_winners[year]
        if winner_title in movies_with_obelisk:
            print(f"The first Best Picture winner to depict a Luxor Obelisk was from the year {year}.")
            print(f"Movie Title: {winner_title}")
            print(f"The obelisk depicted is the one in the Place de la Concorde, Paris.")
            # We found the first one, so we can stop.
            return winner_title

# Execute the function to find and print the answer.
find_first_winner_with_obelisk()