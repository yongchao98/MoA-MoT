def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # Step 1-3: Create a dictionary of battles/campaigns depicted in Lady Butler's
    # paintings and their corresponding years. For events spanning multiple years,
    # the latest year is used for comparison.
    battles = {
        "Battle of Albuera": 1811,
        "Battle of Waterloo": 1815,
        "Battle of Quatre Bras": 1815,
        "Battle of Inkerman": 1854,
        "Battle of Rorke's Drift": 1879,
        "First Boer War": 1881,
        "Battle of Tel el-Kebir": 1882,
        "Nile Expedition": 1885,
    }

    print("Comparing the years of historical events depicted by Lady Butler:")
    # Step 4: Show each event and its year.
    # The loop prints each number in the "equation" for finding the latest year.
    for battle, year in battles.items():
        print(f"Event: {battle}, Year: {year}")

    # Step 5: Find the latest battle using the max() function on the dictionary's values.
    # The 'key=battles.get' argument tells the max function to find the key
    # associated with the maximum value.
    latest_battle_name = max(battles, key=battles.get)
    latest_year = battles[latest_battle_name]

    print("\n--- RESULT ---")
    print(f"The latest historical event depicted in a painting by Lady Butler is the {latest_battle_name}, which took place up to the year {latest_year}.")

find_latest_battle()