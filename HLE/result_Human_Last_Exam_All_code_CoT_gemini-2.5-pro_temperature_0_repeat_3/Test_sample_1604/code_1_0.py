def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the year of the battle depicted.
    # Format: "Painting Title": {"battle": "Battle Name", "year": YYYY}
    paintings = {
        "The Roll Call": {"battle": "Battle of Inkerman (Crimean War)", "year": 1854},
        "Quatre Bras": {"battle": "Battle of Quatre Bras", "year": 1815},
        "Balaclava": {"battle": "Battle of Balaclava (Crimean War)", "year": 1854},
        "The Defence of Rorke's Drift": {"battle": "Battle of Rorke's Drift (Anglo-Zulu War)", "year": 1879},
        "Scotland for Ever!": {"battle": "Battle of Waterloo", "year": 1815},
        "Floreat Etona!": {"battle": "Battle of Laing's Nek (First Boer War)", "year": 1881},
        "The Camel Corps": {"battle": "Battle of Abu Klea (Gordon Relief Expedition)", "year": 1885},
        "A Desperate Stand at Gandamak": {"battle": "Battle of Gandamak (First Anglo-Afghan War)", "year": 1842},
        "The Gate of the Old City": {"battle": "Fall of Baghdad (World War I)", "year": 1917}
    }

    # Find the latest year
    latest_year = 0
    for painting_title in paintings:
        year = paintings[painting_title]["year"]
        if year > latest_year:
            latest_year = year

    # Find the painting and battle for the latest year
    latest_painting_title = ""
    latest_battle_name = ""
    for painting_title, details in paintings.items():
        if details["year"] == latest_year:
            latest_painting_title = painting_title
            latest_battle_name = details["battle"]
            break # Assuming one latest battle for simplicity

    print(f"The latest historical event depicted in a painting by Lady Butler is the {latest_battle_name}.")
    print(f"Painting: '{latest_painting_title}'")
    print(f"Year of Battle: {latest_year}")

find_latest_battle_painting()