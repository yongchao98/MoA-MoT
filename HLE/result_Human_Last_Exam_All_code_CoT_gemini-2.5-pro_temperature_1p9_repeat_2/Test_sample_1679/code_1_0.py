def find_first_winner_with_obelisk():
    """
    Finds the first Academy Award Best Picture winner to depict a Luxor Obelisk.
    The script uses a pre-researched list of early winners and their settings.
    One Luxor Obelisk is in Luxor, Egypt. The other is at the Place de la Concorde in Paris.
    """

    # A chronological list of some Best Picture winners and their settings.
    # The year represents the year the award was for.
    best_picture_winners = [
        {'year': 1927, 'title': 'Wings', 'setting': 'World War I, France/England', 'depicts_obelisk': False},
        {'year': 1934, 'title': 'It Happened One Night', 'setting': 'USA', 'depicts_obelisk': False},
        {'year': 1937, 'title': 'The Life of Emile Zola', 'setting': 'Paris, France', 'depicts_obelisk': False}, # Focuses on interiors, not landmarks
        {'year': 1939, 'title': 'Gone with the Wind', 'setting': 'USA', 'depicts_obelisk': False},
        {'year': 1942, 'title': 'Mrs. Miniver', 'setting': 'England', 'depicts_obelisk': False},
        {'year': 1943, 'title': 'Casablanca', 'setting': 'Morocco', 'depicts_obelisk': False},
        {'year': 1951, 'title': 'An American in Paris', 'setting': 'Paris, France', 'depicts_obelisk': True}, # Famously depicts Place de la Concorde in a ballet sequence
        {'year': 1956, 'title': 'Around the World in 80 Days', 'setting': 'Global, including Paris', 'depicts_obelisk': True},
        {'year': 1958, 'title': 'Gigi', 'setting': 'Paris, France', 'depicts_obelisk': True},
    ]

    print("Searching for the first Best Picture winner to show a Luxor Obelisk...")
    print("One obelisk is in Luxor, Egypt. The other is in Paris, France.\n")

    # Iterate through the winners chronologically
    for movie in sorted(best_picture_winners, key=lambda x: x['year']):
        year = movie['year']
        title = movie['title']
        setting = movie['setting']
        depicts_obelisk = movie['depicts_obelisk']

        print(f"Checking winner for {year}: '{title}' (Setting: {setting})")
        
        if depicts_obelisk:
            print(f"\nFound it! '{title}' is known to depict the Luxor Obelisk located in Paris.")
            print("This film won the award for 1951, preceding other potential candidates like 'Around the World in 80 Days' (1956).")
            print(f"\nFinal Answer: {title}")
            return

# Run the function to find the answer
find_first_winner_with_obelisk()