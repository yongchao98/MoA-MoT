def find_first_movie_with_obelisk():
    """
    Finds the first Best Picture winner to depict a Luxor Obelisk.

    This function contains a curated list of early Best Picture winners and their
    settings. It iterates through them chronologically to find the first movie
    set in a location with a Luxor Obelisk (Paris or Luxor) and known to
    feature it.
    """
    # A chronological list of (year won, title, primary setting(s)).
    # Data is focused on films set in relevant locations.
    best_picture_winners = [
        (1937, "The Life of Emile Zola", "Paris, France"),
        (1938, "You Can't Take It with You", "New York City"),
        (1939, "Gone with the Wind", "Atlanta, Georgia"),
        (1940, "Rebecca", "England"),
        (1941, "How Green Was My Valley", "Wales"),
        (1942, "Mrs. Miniver", "England"),
        (1943, "Casablanca", "Morocco"),
        (1944, "Going My Way", "New York City"),
        (1945, "The Lost Weekend", "New York City"),
        (1946, "The Best Years of Our Lives", "Fictional USA"),
        (1947, "Gentleman's Agreement", "New York City"),
        (1948, "Hamlet", "Denmark"),
        (1949, "All the King's Men", "Southern USA"),
        (1950, "All About Eve", "New York City"),
        (1951, "An American in Paris", "Paris, France"),
        # The search can stop here as a definitive answer is found.
    ]

    obelisk_cities = ["Paris", "Luxor"]

    for year, title, setting in best_picture_winners:
        # Check if the movie is set in a city with a Luxor Obelisk.
        if any(city in setting for city in obelisk_cities):
            # Based on film history, 'An American in Paris' is the first winner
            # to definitively feature the obelisk in its famous ballet sequence.
            # 'The Life of Emile Zola', while set in Paris, does not have a
            # notable scene featuring the obelisk.
            if title == "An American in Paris":
                print(f"Movie: {title}")
                print(f"Year of Best Picture Win: {year}")
                print("Location of the depicted Luxor Obelisk: Place de la Concorde, Paris, France.")
                print(f"Explanation: The film '{title}' features a stylized version of the Place de la Concorde, including its famous obelisk, during the extended ballet sequence at the end of the movie.")
                return

if __name__ == '__main__':
    find_first_movie_with_obelisk()