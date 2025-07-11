import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # A dictionary of Lady Butler's paintings and the historical year of the battle depicted.
    # Format: {"Painting Title": ("Battle Name", Year)}
    paintings = {
        "The 28th Regiment at Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Scotland Forever!": ("Battle of Waterloo", 1815),
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881),
        "The Camel Corps": ("Battle of Abu Klea", 1885),
        "The Gates of Dawn": ("Siege of Mafeking", 1900)
    }

    print("Analyzing Lady Butler's paintings to find the latest depicted battle:")
    print("-" * 70)

    # Print each painting, its depicted battle, and the year.
    for title, (battle, year) in paintings.items():
        print(f"Painting: '{title}' depicts the {battle} which occurred in {year}.")

    # Find the painting with the latest battle year.
    # The `max` function uses a lambda function to specify that the comparison should be based on the year.
    latest_painting_title, (latest_battle_name, latest_year) = max(paintings.items(), key=lambda item: item[1][1])

    print("-" * 70)
    print(f"\nBy comparing the years {', '.join(str(p[1]) for p in paintings.values())}, we can determine the latest event.")
    print(f"\nThe latest historical battle depicted by Lady Butler in this selection is the '{latest_battle_name}', which took place in {latest_year}.")
    print(f"This event was the subject of her painting '{latest_painting_title}'.")

find_latest_battle()