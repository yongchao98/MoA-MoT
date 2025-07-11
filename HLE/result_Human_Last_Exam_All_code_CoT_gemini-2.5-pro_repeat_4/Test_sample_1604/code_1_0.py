import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the year of the battle depicted.
    # Format: "Painting Title": ("Battle Name", Year)
    battle_paintings = {
        "The 28th Regiment at Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Scotland for Ever!": ("Battle of Waterloo", 1815),
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Tel-el-Kebir": ("Battle of Tel El Kebir", 1882),
        "The Charge of the Dorset Yeomanry at Agagia, Egypt": ("Action of Agagia", 1916),
        "Steady the Drums and Fifes": ("Battle of Albuera", 1811),
        "The Fugitives from Lucknow": ("Siege of Lucknow", 1857)
    }

    # Find the painting with the latest battle year
    if not battle_paintings:
        print("No battle paintings found.")
        return

    # The key for max() is a lambda function that accesses the year (index 1) of the dictionary value tuple.
    latest_painting = max(battle_paintings.items(), key=lambda item: item[1][1])
    
    painting_title, (battle_name, battle_year) = latest_painting
    
    print(f"The latest historical battle depicted in a painting by Lady Butler is the '{battle_name}'.")
    print(f"Painting: '{painting_title}'")
    print(f"Battle Year: {battle_year}")

find_latest_battle()