def find_latest_battle():
    """
    Finds the latest historical battle depicted in a selection of Lady Butler's paintings.
    """
    # Step 1, 2 & 3: A dictionary containing painting titles, the battles they depict, and the year of the battle.
    paintings = {
        "The Roll Call": {"battle": "Battle of Inkerman (Crimean War)", "year": 1854},
        "Quatre Bras": {"battle": "Battle of Quatre Bras (Napoleonic Wars)", "year": 1815},
        "The Defence of Rorke's Drift": {"battle": "Battle of Rorke's Drift (Anglo-Zulu War)", "year": 1879},
        "Scotland for Ever!": {"battle": "Battle of Waterloo (Napoleonic Wars)", "year": 1815},
        "Floreat Etona!": {"battle": "Battle of Laing's Nek (First Boer War)", "year": 1881},
        "The Camel Corps": {"battle": "Battle of Abu Klea (Mahdist War)", "year": 1885},
        "Balaclava": {"battle": "Battle of Balaclava (Crimean War)", "year": 1854}
    }

    # Initialize variables to store the latest battle found so far.
    latest_painting_title = None
    latest_battle_info = {"battle": "None", "year": 0}

    print("Analyzing Lady Butler's paintings and their historical battle dates:")
    print("="*60)
    
    # Step 4: Iterate through the data to find the latest year.
    for title, info in paintings.items():
        year = info['year']
        battle = info['battle']
        # This fulfills the requirement to output each number in the "equation"
        print(f"Painting: '{title}', depicts Battle: '{battle}', Year: {year}")
        
        if year > latest_battle_info['year']:
            latest_battle_info = info
            latest_painting_title = title

    # Step 5: Output the final result.
    print("="*60)
    print("\nResult:")
    print(f"The latest historical event depicted among these paintings is the '{latest_battle_info['battle']}'.")
    print(f"This battle, which took place in the year {latest_battle_info['year']}, is depicted in the painting titled '{latest_painting_title}'.")

if __name__ == '__main__':
    find_latest_battle()
    print("\n<<<Battle of Abu Klea (Mahdist War)>>>")