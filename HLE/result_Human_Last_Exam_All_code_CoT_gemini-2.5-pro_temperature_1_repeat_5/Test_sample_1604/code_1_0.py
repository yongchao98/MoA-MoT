def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a known painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the battles they depict.
    # Format: {'Painting Title': ('Battle Name', Year)}
    paintings = {
        'The Roll Call': ('Battle of Inkerman', 1854),
        'Quatre Bras': ('Battle of Quatre Bras', 1815),
        'Balaclava': ('Battle of Balaclava', 1854),
        'The Defence of Rorke\'s Drift': ('Battle of Rorke\'s Drift', 1879),
        'Scotland for Ever!': ('Battle of Waterloo', 1815),
        'Floreat Etona!': ('Battle of Laing\'s Nek', 1881),
        'The Camel Corps': ('Battle of Abu Klea', 1885)
    }

    # Find the painting with the latest battle date
    # We use max() with a lambda function to look at the year (the second item in the tuple value)
    latest_painting_title = max(paintings, key=lambda k: paintings[k][1])
    
    latest_battle_info = paintings[latest_painting_title]
    battle_name = latest_battle_info[0]
    battle_year = latest_battle_info[1]

    print(f"The latest historical battle depicted in a painting by Lady Butler is the '{battle_name}'.")
    print(f"This event took place in the year {battle_year}.")
    print(f"The painting is titled '{latest_painting_title}'.")

find_latest_battle_painting()