def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the year of the battle depicted.
    # Format: "Painting Title": ("Battle Name", year)
    paintings = {
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Scotland Forever!": ("Battle of Waterloo", 1815),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881),
        "A Desperate Stand at Gandamak": ("Battle of Gandamak", 1842),
        "Steady the Drums and Fifes": ("Battle of Albuera", 1811)
    }

    # Find the painting with the latest battle year
    latest_painting = None
    latest_battle_info = (None, float('-inf')) # (Battle Name, Year)

    for painting, (battle, year) in paintings.items():
        if year > latest_battle_info[1]:
            latest_painting = painting
            latest_battle_info = (battle, year)

    battle_name, battle_year = latest_battle_info
    
    print(f"The latest historical battle depicted in a painting by Lady Butler is the {battle_name}.")
    print(f"The painting is titled '{latest_painting}'.")
    print(f"The battle took place in the year {battle_year}.")

find_latest_battle_painting()