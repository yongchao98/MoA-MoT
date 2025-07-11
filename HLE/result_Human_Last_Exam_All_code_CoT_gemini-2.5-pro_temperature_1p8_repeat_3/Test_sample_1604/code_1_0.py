def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler from a predefined list.
    """
    # A dictionary of battles depicted by Lady Butler and the year they occurred.
    # The dictionary keys are the battle names and the values are the years.
    battles = {
        "Battle of Quatre Bras": 1815,
        "Battle of Waterloo": 1815,
        "Battle of Inkerman (depicted in 'The Roll Call')": 1854,
        "Battle of Rorke's Drift": 1879,
        "Battle of Laing's Nek (depicted in 'Floreat Etona!')": 1881,
        "Battle of Tel El Kebir": 1882
    }

    print("Comparing the years of battles depicted by Lady Butler:")
    # Extracting the years into a list for comparison
    years = list(battles.values())
    
    # Printing each battle and its year
    for battle, year in battles.items():
        print(f"- {battle}: {year}")
    
    # Finding the latest year and the corresponding battle
    latest_year = max(years)
    latest_battle = ""
    for battle, year in battles.items():
        if year == latest_year:
            latest_battle = battle
            break # Stop when the first match is found

    # Formatting the numbers for the final equation style output
    years_str = "], [".join(map(str, years))
    
    print("\nTo find the latest battle, we find the maximum of the years:")
    # The "equation" as requested, showing the comparison of all numbers
    print(f"max([{years_str}]) = {latest_year}")

    print("\nThe latest battle depicted is:")
    print(f"{latest_battle}")

find_latest_battle()