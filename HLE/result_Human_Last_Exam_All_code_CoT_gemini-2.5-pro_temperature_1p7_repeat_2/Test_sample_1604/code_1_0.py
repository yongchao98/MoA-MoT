import operator

def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    paintings = {
        "The Remnants of an Army": ("Retreat from Kabul", 1842),
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Scotland for Ever!": ("Battle of Waterloo", 1815),
        "Tel-el-Kebir": ("Battle of Tel El Kebir", 1882),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881)
    }

    # Find the painting with the latest battle year
    latest_painting_title = max(paintings.items(), key=lambda item: item[1][1])
    
    painting_name, (battle_name, battle_year) = latest_painting_title

    print(f"The painting by Lady Butler that depicts the latest historical battle is '{painting_name}'.")
    print(f"It depicts the '{battle_name}', which took place in the year {battle_year}.")

find_latest_battle_painting()