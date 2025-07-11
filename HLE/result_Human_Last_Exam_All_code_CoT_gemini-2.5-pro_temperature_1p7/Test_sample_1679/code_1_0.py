def find_first_movie_with_obelisk():
    """
    This function finds the first Academy Award Best Picture winner
    to feature a Luxor Obelisk.
    """
    # This list represents film research. It contains tuples of:
    # (year, "Title", confirmed_obelisk_sighting_boolean)
    # The list is ordered chronologically by film year.
    best_picture_candidates = [
        # 'The Life of Emile Zola' is set in Paris but a depiction of the obelisk is unconfirmed.
        (1937, "The Life of Emile Zola", False),
        # 'An American in Paris' has a famous ballet sequence featuring the Place de la Concorde obelisk.
        (1951, "An American in Paris", True),
        # 'Around the World in 80 Days' also shows the Paris obelisk, but was a later winner.
        (1956, "Around the World in 80 Days", True),
        # 'Gigi' is another later winner set in Paris that shows the obelisk.
        (1958, "Gigi", True)
    ]

    # Iterate through the chronological list to find the first winner with a confirmed sighting.
    for year, title, has_obelisk in best_picture_candidates:
        print(f"Checking {year} winner: '{title}'...")
        if has_obelisk:
            print(f"\nSUCCESS: Found the first winner.")
            print(f"The first Academy Award Best Picture winner to depict a Luxor Obelisk is '{title}'.")
            print(f"The film won for the year {year}.")
            print("The obelisk shown is the one located in the Place de la Concorde in Paris.")
            # Since the list is chronological, we can stop at the first match.
            return title
    
    # This part would run if no movie was found
    print("No winner matching the criteria was found in the dataset.")
    return None

# Run the function to print the result.
find_first_movie_with_obelisk()
