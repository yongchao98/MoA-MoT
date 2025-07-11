def find_latest_battle():
    """
    This function identifies the latest historical battle depicted in a selection
    of paintings by Lady Butler.
    """
    # A dictionary of Lady Butler's paintings and the years of the battles they depict.
    paintings_and_battles = {
        "The 28th Regiment at Quatre Bras": {"battle": "Battle of Quatre Bras", "year": 1815},
        "Scotland Forever!": {"battle": "Battle of Waterloo", "year": 1815},
        "A Desperate Stand: The Last of the 44th at Gandamak": {"battle": "Last stand at Gandamak", "year": 1842},
        "The Roll Call": {"battle": "Battle of Inkerman (Crimean War)", "year": 1854},
        "The Defence of Rorke's Drift": {"battle": "Battle of Rorke's Drift", "year": 1879},
        "Floreat Etona!": {"battle": "Battle of Laing's Nek", "year": 1881}
    }

    # Find the painting corresponding to the latest battle year.
    latest_painting = None
    latest_year = -1
    years_list = []

    print("Comparing the years of the battles depicted in Lady Butler's paintings:")
    for title, info in paintings_and_battles.items():
        year = info["year"]
        battle = info["battle"]
        years_list.append(year)
        print(f"- '{title}' depicts the {battle}, which occurred in {year}.")
        if year > latest_year:
            latest_year = year
            latest_painting = (title, battle)
    
    # Print the final "equation" showing the comparison of all years.
    print("\nThe equation to find the latest year is determining the maximum value from the set of years:")
    equation = "max("
    for year in sorted(years_list):
        equation += str(year) + ", "
    # Remove the trailing comma and space
    equation = equation[:-2] + ")"
    
    print(f"{equation} = {latest_year}")

    if latest_painting:
        title, battle = latest_painting
        print(f"\nThe latest historical battle depicted is the {battle}, from the year {latest_year}.")
        print(f"This is shown in the painting titled '{title}'.")

if __name__ == "__main__":
    find_latest_battle()