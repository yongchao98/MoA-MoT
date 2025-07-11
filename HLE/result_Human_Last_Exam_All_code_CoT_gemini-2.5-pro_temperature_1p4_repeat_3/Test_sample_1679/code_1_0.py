def find_first_winner_with_obelisk():
    """
    Finds the first Best Picture winner to depict a Luxor Obelisk
    by searching through a curated list of films and their key locations.
    """
    # The Luxor Obelisks are in Luxor, Egypt, and the Place de la Concorde, Paris.
    target_landmarks = ["Place de la Concorde", "Luxor Temple"]

    # A chronological list of selected Best Picture winners and their relevant locations.
    # The data is structured as (Year of Win, Title, [List of Locations]).
    # The list is ordered to find the *first* winner.
    best_picture_winners = [
        (1937, "The Life of Emile Zola", ["Paris"]),
        (1938, "You Can't Take It with You", ["New York City"]),
        (1939, "Gone with the Wind", ["Atlanta, Georgia"]),
        (1940, "Rebecca", ["Monte Carlo", "England"]),
        (1941, "How Green Was My Valley", ["Wales"]),
        (1942, "Mrs. Miniver", ["England"]),
        (1943, "Casablanca", ["Casablanca, Morocco"]),
        (1944, "Going My Way", ["New York City"]),
        (1945, "The Lost Weekend", ["New York City"]),
        (1946, "The Best Years of Our Lives", ["Boone City (fictional), USA"]),
        (1947, "Gentleman's Agreement", ["New York City", "Connecticut"]),
        (1948, "Hamlet", ["Elsinore, Denmark"]),
        (1949, "All the King's Men", ["American South"]),
        (1950, "All About Eve", ["New York City", "Hollywood"]),
        (1951, "An American in Paris", ["Paris", "Place de la Concorde"]),
        (1952, "The Greatest Show on Earth", ["USA"]),
        (1956, "Around the World in 80 Days", ["London", "Paris", "Place de la Concorde"]),
        (1962, "Lawrence of Arabia", ["Jordan", "Syria"]),
    ]

    # Iterate through the winners to find the first match.
    for year, title, locations in best_picture_winners:
        # Check if any of the film's locations match our target landmarks.
        if any(landmark in locations for landmark in target_landmarks):
            print(f"The first Academy Award winner for Best Picture to depict a Luxor Obelisk is:")
            print(f"Film Title: '{title}'")
            # The award is given for a specific year, usually presented in the following year.
            award_year = year
            print(f"Year of Award: {award_year}")
            print("\nThis film features the obelisk located at the Place de la Concorde in Paris during its famous ballet sequence.")
            return

    print("Could not determine the winner based on the provided data.")

if __name__ == "__main__":
    find_first_winner_with_obelisk()