import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a known painting by Lady Butler.
    """
    # A dictionary of battles/events depicted in Lady Butler's paintings and their years.
    # The year for a war is the end year, representing the latest possible time for a depiction.
    battle_paintings = {
        "Battle of Quatre Bras (from 'The 28th Regiment at Quatre Bras')": 1815,
        "Battle of Waterloo (from 'Scotland for Ever!')": 1815,
        "Battle of Inkerman (from 'The Roll Call')": 1854,
        "Battle of Balaclava (from 'Balaclava')": 1854,
        "Battle of Rorke's Drift (from 'The Defence of Rorke's Drift')": 1879,
        "Battle of Laing's Nek (from 'Floreat Etona!')": 1881,
        "Second Boer War (from 'To the Front!')": 1902,
    }

    # Find the battle with the latest year
    # The max function is used on the dictionary's items, comparing them by their value (the year),
    # which is the second element (index 1) of each item tuple.
    latest_battle_item = max(battle_paintings.items(), key=operator.itemgetter(1))

    latest_battle_name = latest_battle_item[0]
    latest_year = latest_battle_item[1]

    # Print the result in a readable format
    print("List of battles/events and their years depicted by Lady Butler:")
    for battle, year in sorted(battle_paintings.items(), key=operator.itemgetter(1)):
        print(f"- {battle}: {year}")

    print("\n--------------------------------------------------")
    print(f"The latest historical event depicted in a painting by Lady Butler is from the:")
    print(f"{latest_battle_name} which concluded in the year {latest_year}.")
    print("--------------------------------------------------")


find_latest_battle()