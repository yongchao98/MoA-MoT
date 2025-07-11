import sys

def find_latest_battle():
    """
    This function identifies the latest battle depicted in a painting by Lady Butler
    from a predefined list.
    """
    # A dictionary mapping the name of the battle to the year it occurred.
    # The names include the painting title for clarity where applicable.
    battles = {
        "Battle of Albuera ('Steady the Drums and Fifes')": 1811,
        "Battle of Waterloo ('Scotland Forever!')": 1815,
        "Battle of Quatre Bras ('Quatre Bras')": 1815,
        "Battle of Balaclava ('Balaclava')": 1854,
        "Battle of Inkerman ('The Roll Call')": 1854,
        "Battle of Rorke's Drift ('The Defence of Rorke's Drift')": 1879,
        "Battle of Tel El Kebir ('Tel-el-Kebir')": 1882,
        "Battle of Abu Klea ('A Desperate Stand at Abu Klea, 1885')": 1885,
    }

    # Find the battle with the maximum year.
    # The max function's key argument specifies a function to be called on each list element
    # prior to making comparisons. Here, it gets the year (the value) for each battle (the key).
    latest_battle_name = max(battles, key=battles.get)
    latest_year = battles[latest_battle_name]

    print("Comparing the years of battles depicted by Lady Butler:")
    # Using a list to ensure consistent output order for the demonstration
    battle_list = sorted(battles.items(), key=lambda item: item[1])
    for battle, year in battle_list:
        print(f"{year}: {battle}")

    print("\n--------------------------------------------------")
    print(f"The latest year is {latest_year}, which corresponds to the {latest_battle_name}.")

    # Writing the final answer to stdout for capture
    # The format requires the answer to be on the last line.
    sys.stdout.write("<<<Battle of Abu Klea>>>")

find_latest_battle()