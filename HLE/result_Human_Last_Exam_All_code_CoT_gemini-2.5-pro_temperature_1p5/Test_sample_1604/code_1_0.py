import operator

def find_latest_battle():
    """
    Finds the latest historical battle among a selection of those painted by Lady Butler.
    """
    # A dictionary of famous paintings/depicted events by Lady Butler and their years.
    # The keys represent the battle or military campaign.
    battles = {
        "Battle of Waterloo (depicted in 'Scotland for Ever!')": 1815,
        "Battle of Quatre Bras (depicted in 'Quatre Bras')": 1815,
        "Battle of Balaclava (depicted in 'Balaclava')": 1854,
        "Crimean War (depicted in 'The Roll Call')": 1855,
        "Battle of Rorke's Drift (depicted in 'The Defence of Rorke's Drift')": 1879,
        "Battle of Laing's Nek (depicted in 'Floreat Etona!')": 1881,
        "Gordon Relief Expedition (depicted in 'A Desert Grave')": 1885,
    }

    print("Analyzing the years of battles painted by Lady Butler:")
    print("-" * 60)

    # Print each battle and its year
    for battle, year in battles.items():
        print(f"Event: {battle}, Year: {year}")

    # Find the battle with the latest year
    if not battles:
        print("The list of battles is empty.")
        return

    # Using max() with a lambda function to find the item with the highest year
    latest_battle, latest_year = max(battles.items(), key=operator.itemgetter(1))

    print("-" * 60)
    print(f"\nResult:")
    print(f"The latest historical event from this list is the {latest_battle}, which occurred in the year {latest_year}.")

find_latest_battle()