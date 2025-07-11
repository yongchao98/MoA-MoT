import collections

def find_latest_battle():
    """
    Identifies the latest historical battle depicted in a painting by Lady Butler
    from a curated list of her works.
    """
    # Step 1: Create a list of tuples containing the painting title, the conflict name, and the year of the event.
    # The year represents the start of the conflict or the specific date of the battle.
    lady_butler_battles = [
        ("The 28th Regiment at Quatre Bras", "Battle of Quatre Bras", 1815),
        ("Scotland for Ever!", "Battle of Waterloo", 1815),
        ("The Roll Call", "Battle of Inkerman", 1854),
        ("Balaclava", "Battle of Balaclava", 1854),
        ("The Defence of Rorke's Drift", "Battle of Rorke's Drift", 1879),
        ("Floreat Etona!", "First Boer War", 1881),
        ("A Detachment of Cavalry in Flanders", "World War I", 1914)
    ]

    print("Analyzing battles depicted by Lady Butler to find the latest one.\n")
    
    # Step 2: Find the entry with the latest year.
    latest_battle_info = max(lady_butler_battles, key=lambda item: item[2])
    
    # Extract the information for the latest battle
    latest_painting_title = latest_battle_info[0]
    latest_battle_name = latest_battle_info[1]
    latest_year = latest_battle_info[2]

    # Step 3: Display the comparison of the years.
    # Get a unique, sorted list of the years for the "equation".
    unique_years = sorted(list(set(item[2] for item in lady_butler_battles)), reverse=True)
    
    print("Comparison of the years of the depicted conflicts:")
    # We "output each number in the final equation" by showing the comparison.
    equation_str = " > ".join(map(str, unique_years))
    print(equation_str)
    
    # Step 4: Print the final answer.
    print(f"\nThe latest historical conflict Lady Butler depicted was {latest_battle_name}, which began in the year {latest_year}.")
    print(f"This is shown in paintings such as '{latest_painting_title}'.")

if __name__ == '__main__':
    find_latest_battle()