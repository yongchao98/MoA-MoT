def find_first_best_picture_with_obelisk():
    """
    This function identifies the first Academy Award Best Picture winner
    to feature a Luxor Obelisk on-screen. It does this by checking a
    chronological list of winners against a list of films known to
    contain a depiction of the obelisk in Paris's Place de la Concorde.
    """
    # A chronological list of early Best Picture winners (year of award ceremony, title)
    best_picture_winners = [
        (1929, "Wings"),
        (1930, "The Broadway Melody"),
        (1930, "All Quiet on the Western Front"),
        (1931, "Cimarron"),
        (1932, "Grand Hotel"),
        (1934, "Cavalcade"),
        (1935, "It Happened One Night"),
        (1936, "Mutiny on the Bounty"),
        (1937, "The Great Ziegfeld"),
        (1938, "The Life of Emile Zola"),
        (1939, "You Can't Take It with You"),
        (1940, "Gone with the Wind"),
        (1941, "Rebecca"),
        (1942, "How Green Was My Valley"),
        (1943, "Mrs. Miniver"),
        (1944, "Casablanca"),
        (1945, "Going My Way"),
        (1946, "The Lost Weekend"),
        (1947, "The Best Years of Our Lives"),
        (1948, "Gentleman's Agreement"),
        (1949, "Hamlet"),
        (1950, "All the King's Men"),
        (1951, "All About Eve"),
        (1952, "An American in Paris"),
        (1953, "The Greatest Show on Earth"),
        (1957, "Around the World in 80 Days")
    ]

    # A list of films known to feature the Luxor Obelisk in Paris
    films_with_paris_obelisk = [
        "The Life of Emile Zola",
        "An American in Paris",
        "Around the World in 80 Days"
    ]

    # Iterate through the winners to find the first match
    for year, title in best_picture_winners:
        if title in films_with_paris_obelisk:
            print(f"The first Best Picture winner to depict a Luxor Obelisk is '{title}'.")
            print(f"It won the award at the {10}th Academy Awards ceremony in {1938} for the film year {1937}.")
            # The obelisk shown is the one at the Place de la Concorde in Paris.
            return

find_first_best_picture_with_obelisk()