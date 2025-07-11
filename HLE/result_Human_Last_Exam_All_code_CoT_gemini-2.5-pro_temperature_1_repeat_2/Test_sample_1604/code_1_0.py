import operator

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # A dictionary of Lady Butler's paintings with their depicted battle and year.
    # Format: {"Painting Title": ("Battle Name", Year)}
    paintings_data = {
        "The 28th Regiment at Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Scotland for Ever!": ("Battle of Waterloo", 1815),
        "The Roll Call": ("Battle of Inkerman", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881),
        "The Camel Corps": ("Battle of Abu Klea", 1885)
    }

    # Find the entry with the latest year
    # The key argument to max() specifies a function to be called on each list element prior to making comparisons.
    # operator.itemgetter(1) gets the value part of the dictionary item.
    # [1] then gets the year from the tuple ("Battle Name", Year).
    latest_painting_entry = max(paintings_data.items(), key=lambda item: item[1][1])

    latest_painting_title = latest_painting_entry[0]
    latest_battle_name, latest_year = latest_painting_entry[1]

    print("Comparing the years of battles depicted by Lady Butler:")
    for title, (battle, year) in sorted(paintings_data.items(), key=lambda item: item[1][1]):
        print(f"- {title}: depicts the {battle} which occurred in {year}.")

    print("\nResult:")
    print(f"The latest historical battle is the '{latest_battle_name}'.")
    print("It was depicted in the painting '{}'.".format(latest_painting_title))
    print(f"The year of the battle is: {latest_year}")

find_latest_battle()