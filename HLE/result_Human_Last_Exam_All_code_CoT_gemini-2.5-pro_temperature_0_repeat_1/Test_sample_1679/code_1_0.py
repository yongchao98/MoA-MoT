def find_first_movie_with_obelisk():
    """
    This function searches through a chronological list of early Academy Award
    Best Picture winners to find the first one known to depict a Luxor Obelisk.
    The two obelisks are in Luxor, Egypt, and Paris, France.
    """
    # Data structure: (Year of Award, "Title", "Primary Setting", Known to Feature Obelisk)
    # The list is ordered chronologically by the award year.
    best_picture_winners = [
        (1928, "Wings", "France/USA", False),
        (1929, "The Broadway Melody", "USA", False),
        (1930, "All Quiet on the Western Front", "Germany", False),
        (1931, "Cimarron", "USA", False),
        (1932, "Grand Hotel", "Germany", False),
        (1933, "Cavalcade", "UK", False),
        (1934, "It Happened One Night", "USA", False),
        (1935, "Mutiny on the Bounty", "Ocean/Tahiti", False),
        (1936, "The Great Ziegfeld", "USA", False),
        (1937, "The Life of Emile Zola", "Paris", False),
        (1938, "You Can't Take It with You", "USA", False),
        (1939, "Gone with the Wind", "USA", False),
        (1940, "Rebecca", "UK", False),
        (1941, "How Green Was My Valley", "UK", False),
        (1942, "Mrs. Miniver", "UK", False),
        (1943, "Casablanca", "Morocco", False),
        (1944, "Going My Way", "USA", False),
        (1945, "The Lost Weekend", "USA", False),
        (1946, "The Best Years of Our Lives", "USA", False),
        (1947, "Gentleman's Agreement", "USA", False),
        (1948, "Hamlet", "Denmark", False),
        (1949, "All the King's Men", "USA", False),
        (1950, "All About Eve", "USA", False),
        (1951, "An American in Paris", "Paris", True),
        (1952, "The Greatest Show on Earth", "USA", False),
        (1953, "From Here to Eternity", "USA", False),
        (1954, "On the Waterfront", "USA", False),
        (1955, "Marty", "USA", False),
        (1956, "Around the World in 80 Days", "Global", True),
    ]

    # Iterate through the list to find the first match
    for year, title, setting, features_obelisk in best_picture_winners:
        if features_obelisk:
            print(f"Searching for the first Best Picture winner featuring a Luxor Obelisk...")
            print(f"Found a match: '{title}' ({year})")
            print(f"This film is set in {setting} and is known to depict the Luxor Obelisk at the Place de la Concorde.")
            return title
    return None

# Run the search and print the final answer
found_movie = find_first_movie_with_obelisk()
if not found_movie:
    print("Could not determine the answer from the provided data.")
