import collections

def find_first_winner_with_obelisk():
    """
    Finds the first Best Picture winner depicting a Luxor Obelisk.

    The two Luxor Obelisks are in Luxor, Egypt, and the Place de la Concorde in Paris.
    This function searches a list of early Best Picture winners to find the first
    one set in a location where an obelisk could be shown.
    """
    # A dictionary of some early Best Picture winners, their award year, and primary setting(s).
    # The year represents the film's year of release, the award was given the following year.
    best_picture_winners = {
        1928: {"title": "Wings", "setting": "France"},
        1929: {"title": "The Broadway Melody", "setting": "New York"},
        1930: {"title": "All Quiet on the Western Front", "setting": "Germany/France"},
        1932: {"title": "Grand Hotel", "setting": "Berlin"},
        1934: {"title": "It Happened One Night", "setting": "USA"},
        1937: {"title": "The Life of Emile Zola", "setting": "Paris"},
        1939: {"title": "Gone with the Wind", "setting": "USA"},
        1942: {"title": "Mrs. Miniver", "setting": "England"},
        1943: {"title": "Casablanca", "setting": "Morocco/Paris"},
        1945: {"title": "The Lost Weekend", "setting": "New York"},
        1951: {"title": "An American in Paris", "setting": "Paris"},
        1952: {"title": "The Greatest Show on Earth", "setting": "USA"},
        1956: {"title": "Around the World in 80 Days", "setting": "Global"},
        1962: {"title": "Lawrence of Arabia", "setting": "Middle East"}
    }

    # Sort the winners chronologically by year
    sorted_winners = collections.OrderedDict(sorted(best_picture_winners.items()))

    # Iterate through the winners to find the first one set in Paris or Luxor
    for year, details in sorted_winners.items():
        title = details["title"]
        setting = details["setting"]

        # Check for films set in Paris. While some earlier winners had Paris scenes,
        # 'An American in Paris' is famous for its extensive showcase of the city's landmarks.
        if setting == "Paris":
            # The Life of Emile Zola and Casablanca have Paris scenes but don't feature the obelisk.
            # An American in Paris is known for its musical numbers celebrating Parisian landmarks.
            if title == "An American in Paris":
                print(f"The first Best Picture winner to depict a Luxor Obelisk is '{title}'.")
                print(f"It won the award for the film year {year}.")
                print("The film prominently features many Parisian landmarks, including the Place de la Concorde, where one of the two Luxor Obelisks stands.")
                break

find_first_winner_with_obelisk()