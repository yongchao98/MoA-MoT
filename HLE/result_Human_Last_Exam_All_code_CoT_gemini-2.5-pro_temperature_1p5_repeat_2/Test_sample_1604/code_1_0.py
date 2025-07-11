def find_latest_battle_painting():
    """
    Analyzes a list of Lady Butler's paintings to find which one
    depicts the latest historical battle.
    """
    # Step 1: Create a list of tuples containing the painting's title,
    # the battle's name, and the year of the battle.
    paintings_data = [
        ("The Roll Call", "Battle of Inkerman", 1854),
        ("The 28th Regiment at Quatre Bras", "Battle of Quatre Bras", 1815),
        ("Scotland Forever!", "Battle of Waterloo", 1815),
        ("The Defence of Rorke's Drift", "Battle of Rorke's Drift", 1879),
        ("The Remnants of an Army", "1842 retreat from Kabul", 1842),
        ("Floreat Etona!", "Battle of Laing's Nek", 1881)
    ]

    # Step 2: Print each painting and its battle year to show the data being analyzed.
    print("Finding the latest battle from the following paintings by Lady Butler:")
    for title, battle, year in paintings_data:
        # This fulfills the requirement to show each number in the 'equation'.
        print(f"- '{title}', depicts event from year: {year}")

    # Step 3: Use the max() function to find the entry with the highest year.
    # The 'key' argument specifies that the comparison should be based on the
    # third element (index 2) of each tuple, which is the year.
    if not paintings_data:
        print("\nNo painting data available to analyze.")
        return

    latest_battle_info = max(paintings_data, key=lambda item: item[2])

    # Step 4: Print the final result.
    final_title, final_battle, final_year = latest_battle_info
    print("\n--- Result ---")
    print(f"The final calculation identifies the latest year as {final_year}.")
    print(f"The latest battle depicted is the '{final_battle}' from the painting '{final_title}'.")


find_latest_battle_painting()
<<<Battle of Laing's Nek>>>