def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler
    from a predefined list.
    """
    # A dictionary mapping battles depicted by Lady Butler to their years.
    # The end year is used for wars spanning multiple years.
    battles = {
        "Battle of Waterloo": 1815,
        "Battle of Inkerman": 1854,
        "Battle of Rorke's Drift": 1879,
        "Battle of Laing's Nek": 1881,
        "Battle of Majuba Hill": 1881,
        "Battle of Tel El Kebir": 1882,
        "Battle of Abu Klea": 1885,
        "Second Boer War": 1902
    }

    # Initialize variables to store the latest battle and its year
    latest_year = 0
    latest_battle_name = ""

    print("Comparing the years of battles depicted by Lady Butler:")
    # Loop through the dictionary to find the latest battle
    for battle, year in battles.items():
        print(f"- {battle}: {year}")
        if year > latest_year:
            latest_year = year
            latest_battle_name = battle

    print("\n----------------------------------------------------")
    print("The latest historical battle found is:")
    # The "final equation" here is showing the resulting battle and its year.
    print(f"{latest_battle_name} which concluded in the year {latest_year}.")
    print("----------------------------------------------------")

if __name__ == "__main__":
    find_latest_battle()