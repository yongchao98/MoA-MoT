import operator

def find_latest_battle():
    """
    This function identifies the latest historical battle depicted in a known
    list of paintings by Elizabeth Thompson, Lady Butler.
    """
    # A dictionary mapping painting titles to a tuple containing the
    # depicted battle/campaign name and its year.
    # For campaigns spanning multiple years, the latest year is used for comparison.
    battles_in_art = {
        "Quatre Bras (1875)": ("Battle of Quatre Bras", 1815),
        "Scotland Forever! (1881)": ("Battle of Waterloo", 1815),
        "The Roll Call (1874)": ("Battle of Inkerman", 1854),
        "Balaclava (1876)": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift (1880)": ("Battle of Rorke's Drift", 1879),
        "Floreat Etona! (1882)": ("Battle of Laing's Nek", 1881),
        "The Camel Corps (1887)": ("Nile Expedition", 1885)
    }

    # Find the entry with the latest year.
    # The max function's key sorts the items by the year, which is the second element
    # in the tuple value (battles_in_art[key][1]).
    latest_entry = max(battles_in_art.items(), key=lambda item: item[1][1])

    latest_painting_title, (latest_battle_name, latest_year) = latest_entry

    print("List of battles depicted by Lady Butler and their years:")
    # Sort by year for clear presentation
    sorted_battles = sorted(battles_in_art.items(), key=lambda item: item[1][1])
    for painting, (battle, year) in sorted_battles:
        print(f"- {battle} ({year}), depicted in '{painting}'")

    print("\n----------------------------------------------------")
    print("Conclusion:")
    print(f"The latest historical military event depicted is the '{latest_battle_name}' which occurred around the year {latest_year}.")
    print(f"It was painted by Lady Butler in her work titled '{latest_painting_title}'.")


if __name__ == '__main__':
    find_latest_battle()