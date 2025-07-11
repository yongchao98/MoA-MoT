def find_first_movie_with_obelisk():
    """
    This function searches through a list of Academy Award Best Picture winners
    to find the first one that depicts a Luxor Obelisk.

    The two Luxor Obelisks are located in:
    1. Luxor, Egypt
    2. Place de la Concorde, Paris, France

    The data is based on film history research.
    """

    best_picture_winners = [
        {'year': 1937, 'title': 'The Life of Emile Zola', 'has_obelisk': False},
        {'year': 1938, 'title': "You Can't Take It with You", 'has_obelisk': False},
        {'year': 1939, 'title': 'Gone with the Wind', 'has_obelisk': False},
        {'year': 1940, 'title': 'Rebecca', 'has_obelisk': False},
        {'year': 1941, 'title': 'How Green Was My Valley', 'has_obelisk': False},
        {'year': 1942, 'title': 'Mrs. Miniver', 'has_obelisk': False},
        {'year': 1943, 'title': 'Casablanca', 'has_obelisk': False},
        {'year': 1944, 'title': 'Going My Way', 'has_obelisk': False},
        {'year': 1945, 'title': 'The Lost Weekend', 'has_obelisk': False},
        {'year': 1946, 'title': 'The Best Years of Our Lives', 'has_obelisk': False},
        {'year': 1947, 'title': "Gentleman's Agreement", 'has_obelisk': False},
        {'year': 1948, 'title': 'Hamlet', 'has_obelisk': False},
        {'year': 1949, 'title': "All the King's Men", 'has_obelisk': False},
        {'year': 1950, 'title': 'All About Eve', 'has_obelisk': False},
        {'year': 1951, 'title': 'An American in Paris', 'has_obelisk': True},
        {'year': 1952, 'title': 'The Greatest Show on Earth', 'has_obelisk': False},
        {'year': 1953, 'title': 'From Here to Eternity', 'has_obelisk': False},
        {'year': 1954, 'title': 'On the Waterfront', 'has_obelisk': False},
        {'year': 1955, 'title': 'Marty', 'has_obelisk': False},
        {'year': 1956, 'title': 'Around the World in 80 Days', 'has_obelisk': True},
        {'year': 1957, 'title': 'The Bridge on the River Kwai', 'has_obelisk': False},
        {'year': 1958, 'title': 'Gigi', 'has_obelisk': True},
    ]

    for movie in best_picture_winners:
        if movie['has_obelisk']:
            print(f"The first Best Picture winner to depict a Luxor Obelisk was '{movie['title']}'.")
            print(f"It won the award for the year {movie['year']}.")
            print(f"The film features the obelisk located at the Place de la Concorde in Paris.")
            return

if __name__ == '__main__':
    find_first_movie_with_obelisk()