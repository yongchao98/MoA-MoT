def find_first_best_picture_with_luxor_obelisk():
    """
    This function identifies the first Best Picture winner
    to feature a Luxor Obelisk and prints the reasoning.
    """
    
    # Step 1: Identify the locations of the two Luxor Obelisks.
    obelisk_location_1 = "Luxor Temple in Luxor, Egypt"
    obelisk_location_2 = "Place de la Concorde in Paris, France"

    # Step 2: Define the film and its winning year.
    film_title = "An American in Paris"
    award_year_for_film = 1951

    # Step 3: Print the reasoning for the answer.
    print("Finding the first Academy Award Best Picture winner to depict a Luxor Obelisk:")
    print("-" * 70)
    print(f"1. There are two Luxor Obelisks: one in {obelisk_location_1} and one in {obelisk_location_2}.")
    print("2. The search of Best Picture winners reveals no films that depict the obelisk in Luxor, Egypt.")
    print("3. The search therefore focuses on the obelisk in Paris.")
    print("4. Early Best Picture winners set in Paris, like 'Casablanca' (1943), do not feature the Place de la Concorde.")
    print(f"5. The film '{film_title}', which won Best Picture for the year {award_year_for_film}, prominently features the Place de la Concorde and its obelisk in the famous final ballet sequence.")
    print("-" * 70)
    
    # Step 4: Print the final answer.
    print(f"The first Best Picture winner to depict a Luxor Obelisk on-screen was: {film_title}")

find_first_best_picture_with_luxor_obelisk()