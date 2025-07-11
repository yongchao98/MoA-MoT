def find_latest_battle():
    """
    Identifies the latest historical battle depicted in a selection of
    Lady Butler's famous paintings.
    """
    # A dictionary mapping painting titles to a tuple containing the
    # battle name and the year the battle occurred.
    paintings_and_battles = {
        "The 28th Regiment at Quatre Bras": ("Battle of Quatre Bras", 1815),
        "Scotland for Ever!": ("Battle of Waterloo", 1815),
        "The Roll Call": ("Battle of Inkerman", 1854),
        "Balaclava": ("Battle of Balaclava", 1854),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Floreat Etona!": ("Battle of Laing's Nek", 1881),
        "After the Battle: Arrival of Lord Wolseley and his Staff at the Bridge of Tel-el-Kebir": ("Battle of Tel El Kebir", 1882)
    }

    print("Comparing the years of battles depicted by Lady Butler:")
    # Loop through the dictionary and print each battle and its year
    for title, (battle, year) in paintings_and_battles.items():
        print(f"- {battle}: {year}")

    # Find the entry with the highest year using a lambda function as a key
    latest_painting_title, (latest_battle, latest_year) = max(
        paintings_and_battles.items(), 
        key=lambda item: item[1][1]
    )

    print("\n-------------------------------------------")
    print("Conclusion:")
    # The final equation is determining the maximum of the years listed above.
    # max(1815, 1815, 1854, 1854, 1879, 1881, 1882) = 1882
    print(f"The latest historical event is the '{latest_battle}' from the year {latest_year}.")
    print(f"This was depicted in her painting titled '{latest_painting_title}'.")

if __name__ == "__main__":
    find_latest_battle()