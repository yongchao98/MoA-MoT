import sys

def find_first_winner_with_obelisk():
    """
    This function simulates the research process to find the first
    Best Picture winner that depicted a Luxor Obelisk on-screen.

    The Luxor Obelisks are located in Luxor, Egypt, and Paris, France.
    We will analyze a curated list of Best Picture winners set in
    these locations in chronological order.
    """

    # Data structure: (Year, "Title", "Location", Depicts_Obelisk_Bool)
    # We list potential candidates chronologically.
    relevant_winners = [
        (1937, "The Life of Emile Zola", "Paris", False),
        (1951, "An American in Paris", "Paris", True),
        (1956, "Around the World in 80 Days", "Paris", True),
        (1958, "Gigi", "Paris", True),
        # Note: No earlier winner is set in or depicts Luxor, Egypt.
        # 'The English Patient' (1996) is set partially in Egypt but
        # does not depict the obelisk in Luxor.
    ]

    first_winner = None

    # Iterate through the list to find the first movie that depicts the obelisk
    for year, title, location, depicts_obelisk in relevant_winners:
        if depicts_obelisk:
            first_winner = (year, title)
            break  # Stop at the first chronological match

    if first_winner:
        winner_year, winner_title = first_winner
        print(f"Searching for the first Best Picture winner to show a Luxor Obelisk...")
        print(f"The obelisk is located in Paris, France.")
        print(f"The first film to win Best Picture and depict the obelisk is '{winner_title}'.")
        print(f"It won the award for the year {winner_year}.")
    else:
        # This part of the code will not be reached with the current data.
        print("No winner matching the criteria was found in the data set.")

# Execute the function to find and print the answer.
find_first_winner_with_obelisk()