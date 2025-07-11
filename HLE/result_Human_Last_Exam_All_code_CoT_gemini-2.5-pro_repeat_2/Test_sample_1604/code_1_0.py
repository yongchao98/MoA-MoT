import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings, the battles they depict, and the year of the battle.
    # Data format: {'Painting Title': ('Battle Name', Year)}
    battle_paintings = {
        'The Remnants of an Army': ('Retreat from Kabul', 1842),
        'The Roll Call': ('Battle of Inkerman', 1854),
        'Balaclava': ('Battle of Balaclava', 1854),
        'Quatre Bras': ('Battle of Quatre Bras', 1815),
        'Scotland Forever!': ('Battle of Waterloo', 1815),
        'The Defence of Rorke\'s Drift': ('Battle of Rorke\'s Drift', 1879),
        'A Desperate Stand': ('Battle of Isandlwana', 1879),
        'Floreat Etona!': ('Battle of Laing\'s Nek', 1881),
        'Tel-el-Kebir': ('Battle of Tel El Kebir', 1882)
    }

    if not battle_paintings:
        print("No battle paintings data found.")
        return

    # Find the painting with the latest battle year
    # We use a lambda function to specify that we want to compare the years (the second element in the tuple value)
    latest_painting_entry = max(battle_paintings.items(), key=lambda item: item[1][1])
    
    painting_title, (battle_name, battle_year) = latest_painting_entry

    print(f"The paintings and their corresponding battle years are:")
    for title, (battle, year) in sorted(battle_paintings.items(), key=lambda item: item[1][1]):
        print(f"- '{title}' depicts the {battle} which occurred in {year}.")

    print("\n--------------------------------------------------")
    print("Conclusion:")
    print(f"The latest historical battle depicted in a painting by Lady Butler is the '{battle_name}'.")
    print(f"It is shown in her painting titled '{painting_title}'.")
    print(f"The battle occurred in the year {battle_year}.")

find_latest_battle()