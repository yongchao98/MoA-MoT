import operator

def find_latest_battle():
    """
    This function finds the latest battle depicted in a selection of
    famous paintings by Lady Butler.
    """
    paintings = [
        {'painting': 'The Roll Call', 'battle': 'Siege of Sevastopol (Crimean War)', 'year': 1855},
        {'painting': 'Quatre Bras', 'battle': 'Battle of Quatre Bras', 'year': 1815},
        {'painting': 'Scotland for Ever!', 'battle': 'Battle of Waterloo', 'year': 1815},
        {'painting': 'The Defence of Rorke\'s Drift', 'battle': 'Battle of Rorke\'s Drift', 'year': 1879},
        {'painting': 'Remnants of an Army', 'battle': 'Retreat from Kabul', 'year': 1842},
        {'painting': 'Balaclava', 'battle': 'Battle of Balaclava', 'year': 1854},
        {'painting': 'Floreat Etona!', 'battle': 'Battle of Laing\'s Nek', 'year': 1881},
        {'painting': 'The Camel Corps', 'battle': 'Battle of Abu Klea', 'year': 1885}
    ]

    print("Comparing the years of battles depicted by Lady Butler:")
    for item in paintings:
        print(f"- {item['battle']}: {item['year']}")
    
    # Find the painting with the latest year
    latest_event = max(paintings, key=operator.itemgetter('year'))

    print("\nResult:")
    print(f"The latest historical battle depicted in a famous painting by Lady Butler is the '{latest_event['battle']}' which occurred in {latest_event['year']}.")
    print(f"The painting is titled '{latest_event['painting']}'.")

if __name__ == "__main__":
    find_latest_battle()