import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of battles depicted by Lady Butler and their corresponding years.
    # The key is the battle/engagement name, and the value is the year it occurred.
    battles = {
        'Battle of Albuera': 1811,
        'Battle of Waterloo': 1815,
        'Battle of Quatre Bras': 1815,
        'Battle of Inkerman': 1854,
        'Battle of Rorke\'s Drift': 1879,
        'Battle of Laing\'s Nek': 1881,
        'Battle of Tel el-Kebir': 1882,
        'Retreat from Mons': 1914,
        'A Zeppelin raid on London': 1915 # The first Zeppelin raid on London was in 1915
    }

    # Find the battle with the latest year.
    # The max function with a key=battles.get finds the key corresponding to the maximum value.
    latest_battle_name = max(battles, key=battles.get)
    latest_year = battles[latest_battle_name]

    print(f"The battles depicted in paintings by Lady Butler include:")
    # Sort the dictionary by year for a chronological output
    sorted_battles = sorted(battles.items(), key=operator.itemgetter(1))
    for battle, year in sorted_battles:
        print(f"- {battle}, which took place in the year {year}")

    print("\nComparing the years...")
    print(f"The latest historical event depicted is the '{latest_battle_name}', which occurred in the year {latest_year}.")

find_latest_battle()