def find_latest_battle():
    """
    Finds the latest historical battle from a list of battles depicted
    in paintings by Lady Butler.
    """
    # A dictionary mapping battle names to their corresponding years.
    # The paintings associated with these battles are by Lady Butler.
    # e.g., 'The Defence of Rorke's Drift' (1879), 'The Camel Corps' (1885), etc.
    battles = {
        "The Battle of Albuera": 1811,
        "The Battle of Quatre Bras": 1815,
        "The Battle of Waterloo": 1815,
        "The Battle of Inkerman": 1854,
        "The Battle of Balaclava": 1854,
        "The Battle of Rorke's Drift": 1879,
        "The Battle of Abu Klea": 1885,
    }

    print("Comparing the years of battles depicted by Lady Butler:")
    # Print each battle and its year for context
    # I am printing each number (year) as requested
    for battle, year in battles.items():
        print(f"- {battle} took place in {year}")
        
    # Find the key (battle name) with the maximum value (year)
    latest_battle_name = max(battles, key=battles.get)
    latest_year = battles[latest_battle_name]

    print("\n-----------------------------------------------------")
    print(f"The latest historical battle depicted in a painting by Lady Butler is:")
    print(f"{latest_battle_name} ({latest_year})")
    print("This battle was portrayed in her painting titled 'The Camel Corps'.")


if __name__ == "__main__":
    find_latest_battle()