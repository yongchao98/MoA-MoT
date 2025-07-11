def find_latest_battle():
    """
    This function identifies the latest historical battle depicted in a selection
    of Lady Butler's paintings.
    """
    # A dictionary mapping painting titles to the year of the event they depict.
    # The key is the painting, and the value is the year.
    battles_by_painting = {
        "The Roll Call": 1854,
        "Quatre Bras": 1815,
        "The Defence of Rorke's Drift": 1879,
        "Scotland for Ever!": 1815,
        "The Remnants of an Army": 1842,
        "Floreat Etona!": 1881,
        "The Camel Corps": 1885
    }

    # A helper dictionary to map painting titles to the specific battle name.
    battle_names = {
        "The Roll Call": "Battle of Inkerman",
        "Quatre Bras": "Battle of Quatre Bras",
        "The Defence of Rorke's Drift": "Battle of Rorke's Drift",
        "Scotland for Ever!": "Battle of Waterloo",
        "The Remnants of an Army": "Retreat from Kabul",
        "Floreat Etona!": "Battle of Laing's Nek",
        "The Camel Corps": "Nile Expedition"
    }

    print("Comparing the years of events depicted in Lady Butler's paintings:")
    all_years = []
    for painting, year in battles_by_painting.items():
        print(f"- '{painting}' depicts an event from {year}")
        all_years.append(year)

    # Find the latest year from the dictionary values
    latest_year = max(battles_by_painting.values())

    # Find the painting(s) associated with that year
    latest_painting = ""
    for painting, year in battles_by_painting.items():
        if year == latest_year:
            latest_painting = painting
            break  # Stop after finding the first match

    latest_battle = battle_names[latest_painting]

    # Display the final calculation as requested
    print("\nFinal Equation:")
    # The 'equation' is finding the maximum value among the years.
    equation_str = f"max({', '.join(map(str, sorted(list(set(all_years)))))}) = {latest_year}"
    print(equation_str)

    print(f"\nThe latest battle depicted is the {latest_battle}, which took place in {latest_year}.")

find_latest_battle()