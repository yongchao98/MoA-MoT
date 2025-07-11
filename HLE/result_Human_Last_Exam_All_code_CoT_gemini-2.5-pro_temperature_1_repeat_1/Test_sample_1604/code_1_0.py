import datetime

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # A list of tuples containing: (Painting Title, Battle/Campaign Name, Year)
    # This list represents some of her most famous works depicting specific, datable events.
    battles_data = [
        ("Steady the Drums and Fifes", "Battle of Albuera", 1811),
        ("The 28th Regiment at Quatre Bras", "Battle of Quatre Bras", 1815),
        ("Scotland for Ever!", "Battle of Waterloo", 1815),
        ("The Roll Call", "Battle of Inkerman", 1854),
        ("Balaclava", "Battle of Balaclava", 1854),
        ("The Defence of Rorke's Drift", "Battle of Rorke's Drift", 1879),
        ("Floreat Etona!", "Battle of Laing's Nek", 1881),
        ("The Camel Corps", "Nile Expedition", 1885),
        ("The Charge of the 21st Lancers at Omdurman", "Battle of Omdurman", 1898)
    ]

    # Find the entry with the maximum year using the max() function.
    # The 'key' argument uses a lambda function to specify that we should look at the third element (the year) of each tuple.
    latest_battle_info = max(battles_data, key=lambda item: item[2])

    painting_title, battle_name, battle_year = latest_battle_info

    print("Analyzing famous paintings by Lady Butler to find the one depicting the latest historical battle.")
    print("--------------------------------------------------")
    print("Considered paintings and their corresponding battle dates:")
    for p, b, y in battles_data:
        print(f"- '{p}' depicts the {b} which took place in {y}")
    print("--------------------------------------------------\n")

    print("The latest historical battle found is:")
    print(f"Painting: '{painting_title}'")
    print(f"Battle: {battle_name}")
    print(f"Year: {battle_year}")

if __name__ == "__main__":
    find_latest_battle()
<<<Battle of Omdurman>>>