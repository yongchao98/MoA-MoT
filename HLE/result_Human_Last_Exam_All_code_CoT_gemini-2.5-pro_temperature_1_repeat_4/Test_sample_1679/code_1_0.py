def find_first_oscar_winner_with_obelisk():
    """
    This function identifies and prints the first Academy Award Best Picture winner
    to feature a Luxor Obelisk.
    """
    # Step 1: There are two Luxor Obelisks. One is in Luxor, Egypt. The other is in the Place de la Concorde in Paris.
    # We must search for the first Best Picture winner featuring either.

    # Step 2: A chronological review of Best Picture winners reveals several films set in Paris.
    # No winners are known to have been filmed at the Luxor Temple in Egypt.
    # Therefore, the search narrows to the obelisk in Paris.

    # Step 3: Potential early winners set in Paris include 'The Life of Emile Zola' (1937)
    # and 'Casablanca' (1942, in flashbacks). However, neither of these films are known
    # for featuring the Place de la Concorde.

    # Step 4: The 1951 Best Picture winner, 'An American in Paris', features a famous 17-minute
    # final ballet sequence. This sequence was designed as a tribute to the city of Paris.
    
    # Step 5: In this sequence, a lavish set was built to represent the Place de la Concorde,
    # with the Luxor Obelisk clearly depicted as the centerpiece. This is the first verifiable
    # instance of the obelisk appearing in a Best Picture-winning film.

    film_title = "An American in Paris"
    award_year = 1951
    location_depicted = "the Place de la Concorde in Paris"

    print("Finding the first Best Picture winner to show a Luxor Obelisk:")
    print(f"The film is '{film_title}'.")
    print(f"It won the Academy Award for Best Picture for the year {award_year}.")
    print(f"It depicted the Luxor Obelisk located in {location_depicted}.")

find_first_oscar_winner_with_obelisk()