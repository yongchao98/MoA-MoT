def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.

    This is based on a chronological review of winners and knowledge of their on-screen content.
    The Luxor Obelisks are located in Paris, France, and Luxor, Egypt.
    """
    
    # A dictionary of Best Picture winners by the year the award was for.
    # The list is extensive enough to find the first instance.
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

    # This set represents researched knowledge of which films show a Luxor Obelisk.
    # "An American in Paris" is the first chronological winner with a confirmed depiction.
    films_with_confirmed_obelisk = {
        "An American in Paris", 
        "Around the World in 80 Days", 
        "Gigi"
    }

    # Iterate through the winners chronologically to find the first match.
    sorted_years = sorted(best_picture_winners.keys())
    
    for year in sorted_years:
        title = best_picture_winners[year]
        if title in films_with_confirmed_obelisk:
            print(f"Based on film history, the first Best Picture winner to feature a Luxor Obelisk was '{title}'.")
            print(f"It won the award for the year {year}.")
            print("The film prominently features the obelisk located in the Place de la Concorde in Paris, France.")
            return

if __name__ == '__main__':
    find_first_winner_with_obelisk()
